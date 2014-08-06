"""Pseudopotential Generators."""
from __future__ import print_function, division

import abc
import os
import tempfile
import collections
import time

from monty.os.path import which
from pseudo_dojo.ppcodes.oncvpsp import OncvOuptputParser
from pymatgen.io.abinitio.pseudos import Pseudo

import logging
logger = logging.getLogger(__name__)


# Possible status of the PseudoGenerator.
_STATUS2STR = collections.OrderedDict([
    (1,  "Initialized"),    # PseudoGenerator has been initialized
    (5,  "Running"),        # PseudoGenerator is running.
    (6,  "Done"),           # Calculation done, This does not imply that results are ok
    (10, "Error"),          #
    (11, "Completed"),      # Execution completed successfully.
])


class Status(int):
    """An integer representing the status of the 'PseudoGenearator`."""
    def __repr__(self):
        return "<%s: %s, at %s>" % (self.__class__.__name__, str(self), id(self))

    def __str__(self):
        """String representation."""
        return _STATUS2STR[self]

    @classmethod
    def as_status(cls, obj):
        """Convert obj into Status."""
        if isinstance(obj, cls):
            return obj
        else:
            # Assume string
            return cls.from_string(obj)

    @classmethod
    def from_string(cls, s):
        """Return a `Status` instance from its string representation."""
        for num, text in _STATUS2STR.items():
            if text == s:
                return cls(num)
        else:
            raise ValueError("Wrong string %s" % s)


class PseudoGenerator(object):
    """
    This object receives a string with the input file and generates a pseudopotential.
    It calls the pp generator in a subprocess to produce results in a temporary directory with a
    non-blocking interface. It also provides an interface to validate/analyze/plot the results 
    produced by the pseudopotential code. Concrete classes must:

        1) call super().__init__() in their constructor.
        2) the object should have the input file stored in self.input_str

    Attributes:
        workdir:
            Working directory (output results are produced in workdir)
        status:
            Flag defining the status of the ps generator.
        retcode:
            Return code of the code
        errors:
            List of strings with errors.
        parser:
            Output parser. None if results are not available because
            the calculations is still running or errors
        results:
            Dictionary with the most important results. None if results are not available because
            the calculations is still running or errors
        pseudo:
            PseudoPotential object. None if not available
    """
    __metaclass__ = abc.ABCMeta

    # Possible status
    S_INIT = Status.from_string("Initialized")
    S_RUN = Status.from_string("Running")
    S_DONE = Status.from_string("Done")
    S_ERROR = Status.from_string("Error")
    S_OK = Status.from_string("Completed")

    ALL_STATUS = [
        S_INIT,
        S_RUN,
        S_DONE,
        S_ERROR,
        S_OK,
    ]

    def __init__(self):
        # Set the initial status.
        self.set_status(self.S_INIT)
        self.errors = []

        # Build a temporary directory
        self.workdir = tempfile.mkdtemp(prefix=self.__class__.__name__)

        # Construct paths for stdin, stdout, stderr and write input file.
        self.stdin_path = os.path.join(self.workdir, "run.in")
        self.stdout_path = os.path.join(self.workdir, "run.out")
        self.stderr_path = os.path.join(self.workdir, "run.err")

    @property
    def status(self):
        return self._status

    @property
    def retcode(self):
        try:
            return self._retcode
        except AttributeError:
            return None

    @property
    def pseudo(self):
        try:
            return self._pseudo
        except AttributeError:
            return None

    @property
    def executable(self):
        return self._executable

    @property
    def input_str(self):
        return self._input_str

    def __repr__(self):
        return "<%s at %s>" % (self.__class__.__name__, os.path.basename(self.workdir))

    def __str__(self):
        return "<%s at %s, status=%s>" % (self.__class__.__name__, os.path.basename(self.workdir), self.status)

    #def __hash__(self):
    #    return hash(self.input_str)
    #def __eq__(self, other):
    #    return self.input_str == other.input_str
    #def __ne__(self, other):
    #    return not self == other

    def start(self):
        """"
        Run the calculation in a sub-process (non-blocking interface)
        Return 1 if calculations started, 0 otherwise.
        """
        if self.status >= self.S_RUN:
            return 0

        logger.info("Running in %s:" % self.workdir)
        with open(self.stdin_path, "w") as fh:
            fh.write(self.input_str)

        # Start the calculation in a subprocess and return.
        args = [self.executable, "<", self.stdin_path, ">", self.stdout_path, "2>", self.stderr_path]
        self.cmd_str = " ".join(args)

        from subprocess import Popen, PIPE
        self.process = Popen(self.cmd_str, shell=True, stdout=PIPE, stderr=PIPE, cwd=self.workdir)
        self.set_status(self.S_RUN, info_msg="Start on %s" % time.asctime)

        return 1

    def poll(self):
        """Check if child process has terminated. Set and return returncode attribute."""
        self._retcode = self.process.poll()

        if self._retcode is not None:
            self.set_status(self.S_DONE)

        return self._retcode

    def wait(self):
        """Wait for child process to terminate. Set and return returncode attribute."""
        self._retcode = self.process.wait()
        self.set_status(self.S_DONE)

        return self._retcode

    def kill(self):
        """Kill the child."""
        self.process.kill()
        self.set_status(self.S_ERROR)
        self.errors.append("Process has beed killed by host code.")
        self._retcode = self.process.returncode

    def set_status(self, status, info_msg=None):
        """
        Set the status.

        Args:
            status:
                Status object or string representation of the status
            info_msg:
                string with human-readable message used in the case of errors (optional)
        """
        assert status in _STATUS2STR

        #changed = True
        #if hasattr(self, "_status"):
        #    changed = (status != self._status)

        self._status = status

        # Add new entry to history only if the status has changed.
        #if changed:
        #    if status == self.S_SUB:
        #        self._submission_time = time.time()
        #        self.history.append("Submitted on %s" % time.asctime())

        #    if status == self.S_OK:
        #        self.history.append("Completed on %s" % time.asctime())

        #    if status == self.S_ABICRITICAL:
        #        self.history.append("Error info:\n %s" % str(info_msg))

        if status == self.S_DONE:
            self.check_status()

        #if status == self.S_OK:
        #    self.on_ok()

        return status

    @abc.abstractmethod
    def check_status(self):
        """
        This function checks the status of the task by inspecting the output and the
        error files produced by the application
        """

    def get_stdin(self):
        return self.input_str

    def get_stdout(self):
        """Returns a string with the stdout of the calculation."""
        if not os.path.exists(self.stdout_path):
            return "Stdout file does not exist"

        with open(self.stdout_path) as out:
            return out.read()

    def get_stderr(self):
        """Returns a string with the stderr of the calculation."""
        if not os.path.exists(self.stdout_path):
            return "Stderr file does not exist"

        with open(self.stderr_path) as err:
            return err.read()

    def rmtree(self):
        """Remove the temporary directory. Return exit status"""
        import shutil
        try:
            shutil.rmtree(self.workdir)
            return 0
        except:
            return 1

    #def on_ok(self):
    #    """
    #    Method called when calculation reaches S_OK
    #    Perform operations to finalize the run. Subclasses should provide their own implementation.
    #    """

    @abc.abstractmethod
    def plot_results(self, **kwargs):
        """Plot the results with matplotlib."""

    def parse_output(self):
        parser = self.OutputParser(self.stdout_path)
        parser.scan()
        return parser

    @property
    def results(self):
        return getattr(self, "_results", None)

    @property
    def plotter(self):
        return getattr(self, "_plotter", None)


class OncvGenerator(PseudoGenerator):
    """
    This object receives an input file for oncvpsp, a string
    that defines the type of calculation (scalar-relativistic, ...)
    runs the code in a temporary directory and provides methods
    to validate/analyze/plot the final results.

    Attributes:
        retcode:
            Retcode of oncvpsp
    """
    from oncvpsp import OncvOuptputParser
    OutputParser = OncvOuptputParser

    def __init__(self, input_str, calc_type):
        super(OncvGenerator, self).__init__()
        self._input_str = input_str
        self.calc_type = calc_type

        calctype2exec = {
            "non-relativistic": which("oncvpspnr.x"),
            "scalar-relativistic": which("oncvpsp.x"),
            "fully-relativistic": which("oncvpspr.x")}

        self._executable = calctype2exec[calc_type]
        if self.executable is None:
            msg = "Cannot find executable for oncvpsp is PATH. Use `export PATH=dir_with_executable:$PATH`"
            raise RuntimeError(msg)

    def check_status(self):
        if self.status == self.S_OK:
            return self._status

        parser = OncvOuptputParser(self.stdout_path)
        try:
            parser.scan()
        except:
            self._status = self.S_ERROR
            return self._status

        logger.info("run_completed:", parser.run_completed)
        if self.status == self.S_DONE and not parser.run_completed:
            logger.info("Run is not completed!")
            self._status = self.S_ERROR

        if parser.run_completed:
            logger.info("setting status to S_OK")
            self._status = self.S_OK
            #########################################
            # Here we initialize results and plotter.
            #########################################
            self._results = parser.get_results()
            self._plotter = parser.make_plotter()

            # Write Abinit pseudopotential.
            filepath = os.path.join(self.workdir, parser.atsym + ".psp8")
            if os.path.exists(filepath): 
                raise RuntimeError("File %s already exists" % filepath)

            # Initialize self.pseudo from file.
            with open(filepath, "w") as fh:
                fh.write(parser.get_pseudo_str())

            self._pseudo = Pseudo.from_file(filepath)
            #print(self.pseudo)

        if parser.ppgen_errors:
            logger.warning("setting status to S_ERROR")
            self._status = self.S_ERROR
            self.errors.extend(parser.ppgen_errors)

        return self._status

    def plot_results(self, **kwargs):
        """Plot the results with matplotlib."""
        if not self.status == self.S_OK:
            logger.warning("Cannot plot results. ppgen status is %s" % self.status)
            return

        # Call the output parser to get the results.
        parser = OncvOuptputParser(self.stdout_path)

        # Build the plotter and plot data according to **kwargs
        plotter = parser.make_plotter()
        plotter.plot_atanlogder_econv()


mock_input = """
# ATOM AND REFERENCE CONFIGURATION
# atsym, z, nc, nv, iexc   psfile
    O    8     1   2   3   psp8
#
# n, l, f  (nc+nv lines)
    1    0    2.0
    2    0    2.0
    2    1    4.0
#
# PSEUDOPOTENTIAL AND OPTIMIZATION
# lmax
    1
#
# l, rc, ep, ncon, nbas, qcut  (lmax+1 lines, l's must be in order)
    0    1.60    0.00    4    7    8.00
    1    1.60    0.00    4    7    8.00
#
# LOCAL POTENTIAL
# lloc, lpopt, rc(5), dvloc0
    4    5    1.4    0.0
#
# VANDERBILT-KLEINMAN-BYLANDER PROJECTORs
# l, nproj, debl  (lmax+1 lines, l's in order)
    0    2    1.50
    1    2    1.00
#
# MODEL CORE CHARGE
# icmod, fcfact
    0    0.0
#
# LOG DERIVATIVE ANALYSIS
# epsh1, epsh2, depsh
   -2.0  2.0  0.02
#
# OUTPUT GRID
# rlmax, drl
    4.0  0.01
#
# TEST CONFIGURATIONS
# ncnf
    3
#
#   nvcnf (repeated ncnf times)
#   n, l, f  (nvcnf lines, repeated follwing nvcnf's ncnf times)
    2
    2    0    2.0
    2    1    3.0
#
    2
    2    0    1.0
    2    1    4.0
#
    2
    2    0    1.0
    2    1    3.0
"""

if __name__ == "__main__":
    pgen = OncvGenerator(input_str=mock_input, calc_type="scalar-relativistic")
    pgen.start()
    pgen.wait()
    print("retcode: %s: " % pgen.retcode)
    #print(pget.get_stdout())
    print(pgen.get_stderr())
    pgen.plot_results()
