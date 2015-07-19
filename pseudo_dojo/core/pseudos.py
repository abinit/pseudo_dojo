# coding: utf-8
"""Public API to access the pseudopotential tables."""
from __future__ import division, print_function, unicode_literals

import os

from pymatgen.core.periodic_table import PeriodicTable
from pymatgen.io.abinitio.pseudos import PseudoTable #, Pseudo


#class DojoTables(object):
#    """This object gathers the different tables in a single namespace."""


class DojoTable(PseudoTable):
    """
    A pseudopotential table provided by pseudo_dojo.
    We subclass `PseudoTable` so that we can easily add 
    extra properties or methods, if needed.
    """
    @classmethod
    def from_dojodir(cls, top):
        """
        Initialize the object for one of the top level directories
        located in pseudo_dojo.pseudos
        """
        # Read metadata from the __init__.py file
        import imp
        module_name = os.path.join(top, "__init__.py")
        meta = imp.load_source(module_name, os.path.join(top, "__init__.py") )

        # Gather all pseudos starting from the current working directory 
        all_symbols = set(element.symbol for element in PeriodicTable().all_elements)
        dirs = [os.path.join(top, d) for d in os.listdir(top) if d in all_symbols]

        paths = []
        for dir in dirs:
            paths.extend(os.path.join(dir, f) for f in os.listdir(dir) 
                         if f.endswith(meta.pseudo_ext))

        return cls(paths).sort_by_z()
