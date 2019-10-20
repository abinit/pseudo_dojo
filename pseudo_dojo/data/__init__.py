"""
Functions providing access to file data for unit tests and tutorials.
Preferred way to import the module is via the import syntax:

import pseudo_dojo.data as pdj_data
"""
import os

from pseudo_dojo.core.pseudos import dojopseudo_from_file, DojoTable


__all__ = [
    "pseudo",
    "pseudos",
    #"ref_file",
    #"ref_files",
]

dirpath = os.path.dirname(__file__)


def pseudopath(filename):
    """Returns the absolute pathname of a pseudo."""
    return os.path.join(dirpath, filename)


def pseudo(filename):
    """Returns a `Pseudo` object."""
    return dojopseudo_from_file(os.path.join(dirpath, filename))


def pseudos(*filenames):
    """Returns a PseudoTable constructed from the input filenames  located in tests/data/pseudos."""
    return DojoTable([dojopseudo_from_file(os.path.join(dirpath, f)) for f in filenames])
