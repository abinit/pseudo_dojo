# coding: utf-8
"""Public API to access the pseudopotential tables."""
from __future__ import division, print_function, unicode_literals

import os
import json

from monty.collections import AttrDict #, dict2namedtuple
from pymatgen.core.periodic_table import PeriodicTable
from pymatgen.io.abinitio.pseudos import PseudoTable 


class DojoTable(PseudoTable):
    """
    A pseudopotential table provided by pseudo_dojo.
    We subclass `PseudoTable` so that we can easily add 
    extra properties or methods, if needed.
    """

    @classmethod
    def from_dojodir(cls, top):
        """
        Initialize the table of pseudos for one of the top level directories
        located in the pseudo_dojo.pseudos directory.

        .. warning::
            
            The table may contain multiple pseudos for a given chemical element.
            Don't use this method unless you need this feature and you know what 
            you are doing.
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

        new = cls(paths).sort_by_z()
        return new

    @classmethod
    def from_djson(cls, djson_file):
        """
        Initialize the table of pseudos for one of **official"" djson files 
        located in the subdirectories on pseudo_dojo.pseudos.

        .. important::

            The table contains one pseudo per element and can be used for production calculations

        djson_file contains the dictionary:

        {
        "dojo_info": {
              "authors": ["M. Giantomassi", "M. J van Setten"],
              "xc_type": "GGA-PBE",
              "pseudo_type": "norm-conserving",
              "generation_date": "Sun Jul 19 21:24:10 CEST 2015",
              "dojo_dir": "ONCVPSP-PBE",
              "description": "String",
              "reference": "paper",
              "tags": ["accuracy", "tag2"]
        },
        "pseudos_metadata": {
            "Si": {
                "basename": "Si-dloc.psp8", 
                "Z_val": 4.0, 
                "l_max": 2, 
                "md5": "ececcf5b26f34676694b630d6bc809ff"
            }, 
            "O": {
                "basename": "O-dmax.psp8", 
                "Z_val": 6.0, 
                "l_max": 2, 
                "md5": "f7d0f3573362d89c81c41fc6b7b3e6ab"
            }
        }, 
        }
        """
        with open(djson_file, "rt") as fh:
            d = json.loads(fh.read())

        dojo_info = AttrDict(**d["dojo_info"])
        #dojo_dir = dojo_info["dojo_dir"]
        meta = d["pseudos_metadata"]

        top = os.path.dirname(djson_file)

        paths, md5dict = [], {}
        for esymb, m in meta.items():
            paths.append(os.path.join(top, esymb, m["basename"]))
            md5dict[m["basename"]] = m["md5"]

        new = cls(paths).sort_by_z()
        new.set_dojo_info(dojo_info)

        errors = new.dojo_check_errors(md5dict)
        if errors:
            raise ValueError("\n".join(errors))

        return new

    @property
    def dojo_info(self):
        try:
            return self._dojo_info
        except AttributeError:
            return {}

    #@dojo_info.setter
    def set_dojo_info(self, dojo_info):
        self._dojo_info = dojo_info

    #def show_dojo_info(self):
    #    print(self.dojo_info)

    def dojo_check_errors(self, md5dict, require_hints=True):
        """
        This function tests whether the table fulfill the requirements 
        imposed by the PseudoDojo. More specifically:

            #. One pseudo per element.

            #. All pseudos should have a valid dojo report with hints

            #. The md5 value computed from the pseudo potential file must agree
               with the one found in the djson file.
        """
        errors = []
        eapp = errors.append

        unique_symbols = set([p.symbol for p in self])
        if len(unique_symbols) != len(self):
            eapp("Found multiple pseudos for a given element.")

        for p in self:
            if not p.has_dojo_report:
                eapp("%s does not have the DOJO_REPORT section" % repr(p))
                continue

            elist = p.dojo_report.check_errors()
            if elist:
                errors.extend(elist)

            if require_hints and not p.dojo_report.has_hints:
                eapp("%s does not have hints" % repr(p))

            if p.md5 != md5dict[p.basename]:
                app("p.mdf5 != mdf5dict[p.basename]\n%s, %s" % (p.md5, md5dict[p.basename]))

        return errors
