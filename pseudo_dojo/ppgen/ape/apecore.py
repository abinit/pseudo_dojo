#!/usr/bin/env python
from __future__ import division, print_function
import os
import sys
import time
import collections
import shutil
from subprocess import Popen, PIPE

import numpy as np

from dftatom.core.qatom import plot_aepp, plot_logders
from dftatom.ape.apeio import (ape_read_waves, ape_read_potentials, ape_read_densities,
                               ape_read_logders, ape_read_dipoles, ape_check_ppeigen, ape_check_ghosts)

__version__ = "0.1"
__status__ = "Development"
__date__ = "$April 26, 2013M$"

##########################################################################################


class ApeError(Exception):
    """Base class for APE Exceptions."""


class ApeWorkDirExists(ApeError):
    """Working directory already exists"""

##########################################################################################

class ApeAeSolver(object):
    """
    This object solves the AE problem using the algorithm and the
    parameters passed via the input generator inpgen.
    """
    _exe = "ape"

    def __init__(self, workdir, inpgen, verbose=0):
        """
        Args:
            workdir:
                path to the working directory.
            inpgen:
                ApeInputGenerator instance.
            verbose:
                Verbosity level.
        """
        self.workdir = os.path.abspath(workdir)
        self.inpgen = inpgen.copy()
        self.verbose = verbose

        self.inpgen.set_calculationmode("ae")

    @classmethod
    def from_template(cls, template, workdir, verbose=0):
        # Create the inpgen to facilitate the modification of the input file.
        inpgen = ApeInputGenerator.from_template(template)
        inpgen.set_calculationmode("ae")

        # Generate new instance
        return cls(workdir, inpgen, verbose=verbose)

    @property
    def exe(self):
        """Executable."""
        return self._exe

    @property
    def input_path(self):
        """Absolute path to the input file."""
        return os.path.join(self.workdir, "ape.in")

    @property
    def output_path(self):
        """Absolute path to the output file."""
        return os.path.join(self.workdir, "ape.out")

    @property
    def stderr_path(self):
        """Absolute path to the stderr file."""
        return os.path.join(self.workdir, "ape.stderr")

    @property
    def aedir(self):
        """Absolute path to the AE directory containing the all-electron results."""
        return os.path.join(self.workdir, "ae")

    # Lazy properties.
    @property
    def ae_waves(self):
        try:
            return self._ae_waves
        except AttributeError:
            self._ae_waves = ape_read_waves(self.aedir)
            return self._ae_waves

    @property
    def ae_potentials(self):
        try:
            return self._ae_potentials
        except AttributeError:
            self._ae_potentials = ape_read_potentials(self.aedir)
            return self._ae_potentials

    @property
    def ae_densities(self):
        try:
            return self._ae_densities
        except AttributeError:
            self._ae_densities = ape_read_densities(self.aedir)
            return self._ae_densities

    def show_input(self, stream=sys.stdout):
        lines  = ["AE INPUT".center(80,"*")]
        lines += self.to_input()
        lines += ["END AE INPUT".center(80,"*")]
        stream.writelines("\n".join(lines)+"\n")

    def to_input(self):
        return self.inpgen.get_strlist()

    def solve(self, **kwargs):
        """Solves the all-electron problem by calling APE in a subprocess."""
        print("Solving the all-electron problem...")
        start = time.time()

        # Write the input file.
        if os.path.exists(self.workdir):
            if kwargs.get("remove_wd", False):
                if self.verbose: print("Will remove %s" % self.workdir)
                shutil.rmtree(self.workdir)
            else:
                raise ApeWorkDirExists(str(self.workdir))

        os.makedirs(self.workdir)
        with open(self.input_path, "w") as fh:
            fh.writelines("\n".join(self.to_input()))

        # Launch APE in a subprocess.
        command = "%s < %s > %s 2> %s" % (self.exe, self.input_path, self.output_path, self.stderr_path)

        process = Popen(command, shell=True, cwd=self.workdir, stderr=PIPE)

        exit_code = process.wait()

        if self.verbose:
            print("AE calculation completed in %.1f s" % (time.time() - start))

        if exit_code != 0:
            with open(self.stderr_path, "r") as stderr:
                raise ApeError("APE returned exit_code %s\n stderr: %s" % (exit_code, stderr.read()))

        with open(self.stderr_path, "r") as stderr:
            error = stderr.read()
            if error:
                raise ApeError("APE exit_code %s\n stderr: %s" % (exit_code, error))

        # TODO Analize the output for possible warnings.
        #self.validate()

        return exit_code

    def plot_waves(self, show=True, savefig=None): 
        plot_aepp(self.ae_waves, show=show, savefig=savefig)

    #def json_dump(self):
    #    json_fname = "__ApeAeSolver__.json"
    #    path = os.path.join(self.workdir, json_fname)
    #    json_pretty_dump(self.to_dict, path)

    #@property
    #def to_dict(self):
    #    "Json-serializable dict representation"
    #    return {
    #        "ae_conf"      : self.ae_conf.to_dict, 
    #        "workdir"      : self.workdir,
    #        "xcfunctional" : self.xcfunctional,
    #        "wave_equation": self.wave_equation,
    #        "radmesh"      : self.radmesh.to_dict,
    #        "control"      : self.control.to_dict,
    #        "@class"       : self.__class__.__name__,
    #        "@module"      : self.__class__.__name__,
    #    }

    #@classmethod
    #def from_dict(cls, d):
    #    "Reconstitute the object from a dict representation  created using to_dict."
    #    #return cls(**{k:v for k,v in d if not k.startswith("@")})

##########################################################################################


class ApePseudoGenerator(object):
    """
    This object uses the data produced by the ApeAeSolver to construct 
    a pseudopotential using the parameters passed to the constructor
    via the input generator inpgen.
    """
    _exe = "ape"

    def __init__(self, workdir, inpgen, ae_solver, verbose=0):
        """
        Args:
            workdir:
                Path to the working directory.
            inpgen:
                ApeInputGenerator instance.
            ae_solver
                ApeAeSolver instance.
            verbose:
                verbosity level.
        """
        self.workdir = os.path.abspath(workdir)
        self.inpgen = inpgen.copy()
        self.ae_solver = ae_solver
        self.verbose = verbose

        self.inpgen.set_calculationmode("pp + pp_test")

    @property
    def exe(self):
        return self._exe

    @property
    def input_path(self):
        return os.path.join(self.workdir, "ape.in")

    @property
    def output_path(self):
        return os.path.join(self.workdir, "ape.out")
                                                        
    @property
    def stderr_path(self):
        return os.path.join(self.workdir, "ape.stderr")

    @property
    def output(self):
        """Output file generated by APE (List of strings)."""
        if not hasattr(self, "_output"):
            with open(self.output_path, "r") as fh:
                self._output = fh.readlines()
        return self._output[:]

    @property
    def ppdir(self):
        """Absolute path to the ae directory containing the all-electron results."""
        return os.path.join(self.workdir, "pp")

    @property
    def testsdir(self):
        """Absolute path to the directory containing test results."""
        return os.path.join(self.workdir, "tests")

    @property
    def ae_waves(self):
        return self.ae_solver.ae_waves

    # Lazy properties.
    @property
    def pp_waves(self):
        try:
            return self._pp_waves
        except AttributeError:
            self._pp_waves = ape_read_waves(self.ppdir)
            return self._pp_waves
                                                       
    @property
    def pp_potentials(self):
        try:
            return self._pp_potentials
        except AttributeError:
            self._pp_potentials = ape_read_potentials(self.ppdir)
            return self._pp_potentials
                                                       
    @property
    def pp_densities(self):
        try:
            return self._pp_densities
        except AttributeError:
            self._pp_densities = ape_read_densities(self.ppdir)
            return self._pp_densities

    @property
    def ae_logders(self):
        try:
            return self._ae_logders
        except AttributeError:
            self._ae_logders, self._pp_logders = ape_read_logders(self.testsdir)
            return self._ae_logders

    @property
    def pp_logders(self):
        try:
            return self._pp_logders
        except AttributeError:
            self._ae_logders, self._pp_logders = ape_read_logders(self.testsdir)
            return self._pp_logders

    @property
    def ghosts(self):
        """Ghost states. Evaluates to True if ghost states are found."""
        try:
            return self._ghosts
        except AttributeError:
            self._ghosts = self._check_ghosts()
            return self._ghosts

    @property
    def dipoles(self):
        """Dipole matrix elements."""
        try:
            return self._dipoles
        except AttributeError:
            self._dipoles = ape_read_dipoles(self.output)
            return self._dipoles

    def path_in_workdir(self, filename):
        """Generates the absolute path in the working directory."""
        return os.path.join(self.workdir, filename)

    def show_input(self, stream=sys.stdout):
        lines  = ["PP INPUT".center(80,"*")]
        lines += self.to_input()
        lines += ["END PP INPUT".center(80,"*")]
        stream.writelines("\n".join(lines)+"\n")

    def to_input(self):
        """Returns the APE input (list of strings)"""
        return self.inpgen.get_strlist()

    def pseudize(self, **kwargs):
        """Performs the pseudization by calling APE in a subprocess."""
        print("Starting pseudization...")
        start = time.time()

        # Write the input file.
        if os.path.exists(self.workdir):
            if kwargs.get("remove_wd", False):
                if self.verbose: print("Will remove %s" % self.workdir)
                shutil.rmtree(self.workdir)
            else:
                raise ApeWorkDirExists(str(self.workdir))
        os.makedirs(self.workdir)

        # FIXME: Have to copy data dir
        #shutil.copy(os.path.join(self.ae_solver.aedir, "data"), self.workdir)
        shutil.copytree(self.ae_solver.aedir, os.path.join(self.workdir, "ae"))

        with open(self.input_path, "w") as fh:
            fh.writelines("\n".join(self.to_input()))

        # Launch APE in a subprocess.
        command = "%s < %s > %s 2> %s" % (self.exe, self.input_path, self.output_path, self.stderr_path)

        process = Popen(command, shell=True, cwd=self.workdir, stderr=PIPE)

        exit_code = process.wait()

        if self.verbose:
            print("Pseudization completed in %.1f s" % (time.time() - start))

        if exit_code != 0:
            with open(self.stderr_path, "r") as stderr:
                raise ApeError("APE exit_code %s\n stderr: %s" % (exit_code, stderr.read()))

        with open(self.stderr_path, "r") as stderr:
            error = stderr.read()
            if error:
                raise ApeError("APE exit_code %s\n stderr: %s" % (exit_code, error))

        # TODO Analize the output for possible warnings.

        # Check ghost-states 
        self._check_ghosts()

        # Check PP eigenvalues
        self._check_ppeigen()                                                

        # Check logarithmic derivative.

        # Dump quality control file.
        #with open(self.path_in_workdir("quality.json"), "w") as fh:
        #    json.dump(self.quality_control(), fh, indent=4, sort_keys=4)

        return exit_code

    #def jsone_dump(self, protocol=0):
    #    path = os.path.join(self.workdir, self.json_fname)
    #    with open(path, "w") as fh:
    #        json.dump(self, fh, protocol=protocol)

    def _check_ppeigen(self):
        return ape_check_ppeigen(self.output)

    def _check_ghosts(self):
        """Check for presence of ghost states. Reads data from ape.out."""
        return ape_check_ghosts(self.output)

    def check_logders(self):
        """Check the quality of the log derivatives."""
        merits = {}
        for (state, pp_ld) in self.pp_logders.items():
            ae_ld = self.ae_logders[state]
            rmesh = ae_ld.rmesh
            f = np.abs(np.tan(ae_ld.values) - np.tan(pp_ld.values))
            from scipy.interpolate import UnivariateSpline
            spline = UnivariateSpline(rmesh, f, s=0)
            merits[state] = spline.integral(rmesh[0], rmesh[-1])  / (rmesh[-1] - rmesh[0])

        return merits

    #def quality_control(self):
    #    return {
    #        "ghosts"        : self.ghosts,
    #        "ppeig_quality" : self._check_ppeigen,
    #        "logder_quality": self.check_logders(),
    #    }

    def plot_waves(self, show=True, savefig=None): 
        plot_aepp(self.ae_waves, pp_funcs=self.pp_waves, show=show, savefig=savefig)

    def plot_logders(self, show=True, savefig=None): 
        plot_logders(self.ae_logders, self.pp_logders, show=show, savefig=savefig)

    #def plot_potentials(self, vname):
    #    plot_potentials(self.ae_pots, pp_pots=self.pp_pots, show=show, savefig=savefig)

    #def plot_densities(self):
    #    plot_densities(self.ae_densities, pp_dens=self.pp_densities, show=show, savefig=savefig)

    #def plot_core_correction(self):

##########################################################################################


def ape_parse_input(input): 
    """
    Parses an APE input file, returns an ordered dictionary {varname: varvalue}

    Args:
        input:
            filename of list of lines.
    
    .. note: varnames are all lower case, the name of variables corresponding
             to APE blocks starts with the character %
    """
    if isinstance(input, str):
        with open(input, "r") as fh:
            input = fh.readlines()

    ovars = collections.OrderedDict()
                                                                         
    # Remove comments and empty lines.
    lines =[]
    for line in input:
        l = line.strip() 
        if l and not l.startswith("#"): 
            lines.append(l)
                                                                         
    # Parse the file: treat blocks and scalar variables separately.
    while True:
        try:
            line = lines.pop(0)
        except IndexError:
            break
                                                                         
        if line.startswith("%"):
            # Not that varname contains the char "%" so that
            # we can separate scalar variables from blocks.
            # print("Begin block: %s" % line)
            var_name = line.lower().strip()
            assert var_name not in ovars
            ovars[var_name] = []
                                                                         
            # Consume the block.
            while True:
                l = lines.pop(0)
                if l.startswith("%"):
                    break
                else:
                    ovars[var_name].append(l)
        else:
            #print("Got scalar: %s" % line)
            var_name, var_value = line.split("=")
            var_name = var_name.lower().strip()
            assert var_name not in ovars
            ovars[var_name] = var_value
                                                                         
    return ovars


class ApeInputGenerator(object):
    """
    This object facilitates the creation/modification of APE inputs.
    An ApeInputGenerator has a template (list of strings) and new variables (list of strings)
    that can be added at run-time via the set methods.

    The method get_strlist returns a new APE input where the new variables replace the ones
    defined in the template. The preferred way to generate an instance is via the class method from_template.
    """
    def __init__(self, template=None, newvars=None):
        self._template = template.copy() if template else {}
        self._newvars = newvars.copy() if newvars else {}

    def __str__(self):
        return "\n".join(self.get_strlist())

    @classmethod
    def from_template(cls, template):
        """Initialize the object from a template (string or list of strings)."""
        if isinstance(template, str):
            with open(os.path.abspath(template), "r") as fh:
                template = ape_parse_input(fh.readlines())
        return cls(template=template)

    #@classmethod
    #def from_objects(cls, *objs)
    #    template = collections.OrderedDict()
    #    for o in objs:
    #        template.update(o.to_vars())
    #    return cls(template=template)

    @property
    def varnames(self):
        """Returns a set with the variable names."""
        return set(list(self._newvars) + list(self._template))

    def copy(self):
        """Shallow copy of self."""
        return ApeInputGenerator(template=self._template, newvars=self._newvars)

    def reset(self):
        """Reset the inpgen by emptying the list of new variables."""
        self._newvars = collections.OrderedDict()

    def _is_scalar(self, varname):
        """True if varname is a scalar variable."""
        return not varname.startswith("%")

    def get_scalar(self, varname):
        """Returns the value of a scalar variable."""
        varname = varname.lower()
        try:
            return self._newvars[varname]
        except KeyError:
            return self._template[varname]

    def set_scalar(self, varname, value):
        """Set the (new) value of a scalar variable"""
        self._newvars[varname.lower()] = value

    def get_block(self, varname):
        """Returns the value of a block variable."""
        if varname[0] != "%": varname = "%" + varname
        varname = varname.lower()
        try:
            return self._newvars[varname]
        except KeyError:
            return self._template[varname]

    def set_block(self, varname, value):
        """Set the (new) value of a block variable."""
        if varname[0] != "%": varname = "%" + varname
        self._newvars[varname.lower()] = value

    def set_calculationmode(self, mode):
        """Change the calculation mode."""
        self.set_scalar("calculationmode", mode)

    def set_llocal(self, llocal):
        """Change the value of Llocal (l for the local part)."""
        self.set_scalar("llocal", llocal)

    def set_correction(self, corecorrection):
        """Change the value of CoreCorrection"""
        self.set_scalar("corecorrection", corecorrection)

    #def set_ppcomponents(self, *strings)

    def get_strlist(self):
        """Return a list of strings with the APE input."""
        # Copy the template and update the keys with those in _newvars
        d = self._template.copy()
        d.update(self._newvars)

        # Format string depending of the type (scalar, block).
        lines = []
        for (k,v) in d.items():
            if self._is_scalar(k):
                lines += ["%s = %s " % (k, v)]
            else:
                lines += ["%s" % k]
                for l in v:
                    lines += ["%s" % l]
                lines += ["%"]
        return lines

    #@abc.abstractproperty
    #def Z(self):
    #    "The atomic number of the atom."

    #@abc.abstractproperty
    #def Z_val(self):
    #    "Valence charge"

    #@property
    #def element(self):
    #    "Pymatgen Element"
    #    return _periodic_table[self.Z]

    #@property
    #def symbol(self):
    #    "Element symbol."
    #    return self.element.symbol

    #@abc.abstractproperty
    #def l_max(self):
    #    "Maximum angular momentum"

    #@abc.abstractproperty
    #def l_local(self):
    #    "Angular momentum used for the local part"

    #@abc.abstractproperty
    #def xc_family(self):
    #    "XC family e.g LDA, GGA, MGGA."

    #@abc.abstractproperty
    #def xc_flavor(self):
    #    "XC flavor e.g PW, PW91, PBE."

    #@property
    #def xc_type(self):
    #    "XC identifier e.g LDA-PW, GGA-PBE, GGA-revPBE"
    #    return "-".join([self.xc_family, self.xc_flavor])

    #@property
    #def isnc(self):
    #    "True if norm-conserving pseudopotential"
    #    return isinstance(self, NcPseudo)

    #@property
    #def ispaw(self):
    #    "True if PAW pseudopotential"
    #    return isinstance(self, PawPseudo)

    #@abc.abstractproperty
    #def has_soc(self):
    #    "True if pseudo contains spin-orbit coupling."

    #@abc.abstractmethod
    #def num_of_projectors(self, l='s'):
    #    "Number of projectors for the angular channel l"

    #@abc.abstractmethod
    #def generation_mode
    #    "scalar scalar-relativistic, relativistic"

    #@abc.abstractproperty
    #def nlcc_radius(self):
    #    """
    #    Radius at which the core charge vanish (i.e. cut-off in a.u.).
    #    Returns 0.0 if nlcc is not used.
    #    """

    #@property
    #def has_nlcc(self):
    #    "True if the pseudo is generated with non-linear core correction."
    #    return self.nlcc_radius > 0.0

    #@abc.abstractproperty
    #def paw_radius(self):
    #    "Radius of the PAW sphere in a.u."
