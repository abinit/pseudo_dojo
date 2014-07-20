"""Pseudopotential Generators."""
from __future__ import print_function, division

from monty.os.path import which
from .oncvpsp import OncvOuptputParser

import abc


class PseudoGenerator(object):
    """
    This object receives a string with the input file and generates a pseudopotential.
    It calls the pp generator in a subprocess to produce results in a temporary directory with a
    non-blocking interface implemented here.
    It also provides an interface to validate/analyze/plot the final results produced by the code.
    The implementation of the constructor is delegated to the sub-class, the only requirement is
    that the object should have the input file stored in self.input_str

    Attributes:
        retcode:
            Return code of the code
    """
    __metaclass__ = abc.ABCMeta

    @property
    def retcode(self):
        try:
            return self._retcode
        except AttributeError:
            return None

    @property
    def executable(self):
        return self._executable

    @property
    def input_str(self):
        return self._input_str

    def start(self):
        """"
        Run the calculation in a sub-process (non-blocking interface)
        """
        # Build a temporary directory
        import tempfile
        import os
        self.dirpath = tempfile.mkdtemp(prefix="ONVCV_RUN_")
        print("Running in %s:" % self.dirpath)

        # Construct paths for stdin, stdout, stderr and write input file.
        self.stdin_path = os.path.join(self.dirpath, "run.in")
        self.stdout_path = os.path.join(self.dirpath, "run.out")
        self.stderr_path = os.path.join(self.dirpath, "run.err")

        with open(self.stdin_path, "w") as fh:
            fh.write(self.input_str)

        # Start the calculation in a subprocess and return.
        args = [self.executable, "<", self.stdin_path, ">", self.stdout_path, "2>", self.stderr_path]
        self.cmd_str = " ".join(args)

        from subprocess import Popen, PIPE
        self.process = Popen(self.cmd_str, shell=True, stdout=PIPE, stderr=PIPE, cwd=self.dirpath)

    def poll(self):
        """Check if child process has terminated. Set and return returncode attribute."""
        self._retcode = self.process.poll()
        return self._retcode

    def wait(self):
        """Wait for child process to terminate. Set and return returncode attribute."""
        self._retcode = self.process.wait()
        return self._retcode

    def kill(self):
        """Kill the child."""
        self.process.kill()
        self._retcode = self.process.returncode

    def get_stdout(self):
        """Returns a string with the stdout of the calculation."""
        with open(self.stdout_path) as out:
            return out.read()

    def rmtree(self):
        """Remove the temporary directory. Return exit status"""
        import shutil
        try:
            shutil.rmtree(self.dirpath)
            return 0
        except:
            return 1

    def get_stderr(self):
        """Returns a string with the stderr of the calculation."""
        with open(self.stderr_path) as err:
            return err.read()

    @abc.abstractmethod
    def plot_results(self, **kwargs):
        """Plot the results with matplotlib."""


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
    def __init__(self, input_str, calc_type):
        self._input_str = input_str
        self.calc_type = calc_type

        calctype2exec = {
            "non-relativistic": which("oncvpspnr.x"),
            "scalar-relativistic": which("oncvpsp.x"),
            "fully-relativistic": which("oncvpspr.x")}

        self._executable = calctype2exec[calc_type]
        if self.executable is None:
            msg = "Cannot find executable %s is PATH\n Use export PATH=/dir_with_exec:$PATH" % self.executable
            raise RuntimeError(msg)

    def plot_results(self, **kwargs):
        """Plot the results with matplotlib."""
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
    pgen = OncvGenerator(mock_input, "scalar-relativistic")
    pgen.start()
    pgen.wait()
    print("retcode: %s: " % pgen.retcode)
    #print(pget.get_stdout())
    print(pgen.get_stderr())
    pgen.plot_results()