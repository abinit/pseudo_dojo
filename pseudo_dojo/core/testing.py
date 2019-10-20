# coding: utf-8
# flake8: noqa
"""
Common test support for pseudo_dojo.

This single module should provide all the common functionality for tests
in a single location, so that test scripts can just import it and work right away.
"""
from monty.os.path import which
from pymatgen.util.testing import PymatgenTest

import numpy.testing.utils as nptu
import subprocess
import tempfile


def cmp_version(this, other, op=">="):
    """
    Compare two version strings with the given operator `op`

    >>> assert cmp_version("1.1.1", "1.1.0") and not cmp_version("1.1.1", "1.1.0", op="==")
    """
    from pkg_resources import parse_version
    from monty.operator import operator_from_str
    op = operator_from_str(op)
    return op(parse_version(this), parse_version(other))


def has_abinit(version=None, op=">="):
    """
    True if abinit is in $PATH.
    If version is not None, abinit version op version is evaluated and the result is returned.
    False if condition is not fulfilled or the execution of `abinit -v` raised CalledProcessError
    """
    abinit = which("abinit")
    if abinit is None: return False
    if version is None: return abinit is not None

    try:
        abinit_version = str(subprocess.check_output(["abinit", "-v"]))

    except subprocess.CalledProcessError:
        # Some MPI implementations require the mpirunner.
        try:
            abinit_version = subprocess.check_output(["mpirun", "-n", "1", "abinit", "-v"])
        except subprocess.CalledProcessError:
            try:
                abinit_version = subprocess.check_output(["mpiexec", "-n", "1", "abinit", "-v"])
            except subprocess.CalledProcessError as exc:
                #logger.warning(exc.output)
                return False

    return cmp_version(abinit_version, version, op=op)


def has_matplotlib(version=None, op=">="):
    """
    True if matplotlib is installed.
    If version is None, the result of matplotlib.__version__ `op` version is returned.
    """
    try:
        #have_display = "DISPLAY" in os.environ
        import matplotlib
        matplotlib.use("Agg")  # Use non-graphical display backend during test.

    except ImportError:
        print("Skipping matplotlib test")
        return False

    # http://stackoverflow.com/questions/21884271/warning-about-too-many-open-figures
    import matplotlib.pyplot as plt
    plt.close("all")

    if version is None: return True

    return cmp_version(matplotlib.__version__, version, op=op)


def has_seaborn():
    """True if seaborn is installed."""
    try:
        import seaborn.apionly as sns
        return True
    except ImportError:
        return False


class PseudoDojoTest(PymatgenTest):
    """Extends PymatgenTest with PseudoDojo-specific methods """

    @staticmethod
    def assert_almost_equal(actual, desired, decimal=7, err_msg='',
                            verbose=True):
        """
        Alternative naming for assertArrayAlmostEqual.
        """
        return nptu.assert_almost_equal(actual, desired, decimal, err_msg, verbose)

    @staticmethod
    def assert_equal(actual, desired, err_msg='', verbose=True):
        """
        Alternative naming for assertArrayEqual.
        """
        return nptu.assert_equal(actual, desired, err_msg=err_msg, verbose=verbose)

    @staticmethod
    def which(program):
        """Returns full path to a executable. None if not found or not executable."""
        return which(program)

    @staticmethod
    def has_abinit(version=None, op=">="):
        """Return True if abinit is in $PATH and version is op min_version."""
        return has_abinit(version=None, op=op)

    @staticmethod
    def has_matplotlib(version=None, op=">="):
        return has_matplotlib(version=version, op=op)

    @staticmethod
    def has_seaborn():
        """True if seaborn is installed."""
        return has_seaborn()

    def get_tmpname(self, **kwargs):
        """Invoke mkstep with kwargs, return the name of a temporary file."""
        fd, tmpname = tempfile.mkstemp(**kwargs)
        return tmpname

    def has_nbformat(self):
        """Return True if nbformat is available and we can test the generation of ipython notebooks."""
        try:
            import nbformat
            return True
        except ImportError:
            return False
