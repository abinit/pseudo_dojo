from __future__ import division, print_function, unicode_literals

import os

from collections import defaultdict, OrderedDict, Mapping
from monty.dev import deprecated
from monty.functools import lazy_property
from monty.design_patterns import singleton
from pseudo_dojo.core.pseudos import DojoTable
from pseudo_dojo.pseudos import dojotable_absdir


@deprecated(message="Official PseudoDojo API will be released soon!")
def get_pseudos(top):
    """
    Find pseudos within top, return :class:`PseudoTable` object sorted by atomic number Z.
    """
    from monty.os.path import find_exts
    from pymatgen.io.abinit.pseudos import PseudoTable, Pseudo
    exts = ("psp8",)
    pseudos = []
    for p in find_exts(top, exts, exclude_dirs="_*"):
        try:
            pseudos.append(Pseudo.from_file(p))
        except Exception as exc:
            from warnings import warn
            warn("Exception in pseudo %s:\n%s" % (p.filepath, exc))

    return PseudoTable(pseudos).sort_by_z()


class TableMetadata(object):
    """
    Metadata associated to one of the official `PseudoDojo` tables.
    """
    def __init__(self, table_dir, djson_name):
        """
        Args:
            table_dir: basename of the directory containing the pseudos
            djson_name: name of the json file in `table_dir` with the
                list of pseudos and metatada.
        """
        self.table_dir = table_dir
        self.dojo_absdir = dojotable_absdir(table_dir)
        self.djson_name = djson_name
        self.djson_path = os.path.join(self.dojo_absdir, djson_name)


@singleton
class OfficialTables(Mapping):
    """
    Official tables provided by the PseudoDojo project.
    """
    def __init__(self):
        self._tables = tables = OrderedDict()

        # Naming scheme: [PP_TYPE]-[XC_NAME]-[TABLE_NAME_WITH_VERSION]-[DJSON_NAME]
        #
        # where:
        #
        #  PP_TYPE in ("PAW", "NC") defines the pseudos type.
        #  XC_NAME gives the functional type e.g. PBE
        #  TABLE_NAME_WITH_VERSION is the name of the table e.g. ONCVPSP, JTH ...
        #  DJSON_NAME is the name of the djson file located in the pseudodojo directory.
        #  In the case of NC pseudos, the SO term is treated via djson e.g. "accuracy-so.djson"
        #  maybe FR|SR for fully or scalar-relativistic?
        #
        tables["ONCVPSP-PBE-PDv0.2-accuracy"] = TableMetadata("ONCVPSP-PBE-PDv0.2", "accuracy.djson")

        #tables["ONCV-GYGYv0.2-PBE"] = TableMetadata("ONCVPSP-PBE", "htc.djson")
        #tables["PAW-JTHv0.2-PBE"] = TableMetadata("ONCVPSP-PBE", "accuracy.djson")
        #tables["PAW-GBRVv0.2-PBE"] = TableMetadata("ONCVPSP-PBE", "efficiency.djson")

    # ABC protocol.
    def __iter__(self):
        return self._tables.__iter__()

    def __len__(self):
        return self._tables.__len__()

    def __getitem__(self, key):
        """Called by self[key]"""
        v = self._tables[key]

        # Return v if v is table else parse the files, construct table and store it.
        if not isinstance(v, TableMetadata):
            return v

        new_table = DojoTable.from_djson_file(v.djson_path)
        new_table.dojo_name = key
        self._tables[key] = new_table

        return new_table

    def select_tables(self, pp_type=None, xc=None):
        """
        """
        tables = self.values()
        tables = [table for table in tables if table.dojo_info.pp_type == pp_type]
        tables = [table for table in tables if table.dojo_info.xc == xc]
        return tables

    #@lazy_property
    #def available_xcs(self):
    #    return set(t.xc for t in self.values())

    def all_nctables(self, xc=None):
        """Return the list of norm-conserving tables."""
        tables = [table for table in self.values() if table.dojo_info.isnc]
        if xc is not None:
            tables = [t for t in tables if t.xc == xc]
        return tables

    def all_pawtables(self, xc=None):
        """Return the list of PAW tables."""
        tables = [table for table in self.values() if table.dojo_info.ispaw]
        if xc is not None:
            tables = [t for t in tables if t.xc == xc]
        return tables

    def select_pseudos(self, symbol, pp_type=None, xc=None):
        """
        Return the full list of Pseudo objects available in the DojoTables 
        with the given `pp_type` and XC functional `xc`.
        """
        tables = self.select_tables(pp_type=pp_type, xc=xc)
        pseudos = []
        for tab in tables:
            pseudos.extend(tab.pseudo_with_symbol(symbol, allow_multi=True))
        return pseudos

    def select_ncpseudos(self, symbol, xc=None):
        return self.select_pseudos(symbol, pp_type="NC", xc=xc)

    def select_pawpseudos(self, symbol, xc=None):
        return self.select_pseudos(symbol, pp_type="PAW", xc=xc)
