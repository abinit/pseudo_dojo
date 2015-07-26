from __future__ import division, print_function, unicode_literals

import os

from monty.dev import deprecated


@deprecated(message="Official PseudoDojo API will be released soon!")
def get_pseudos(top):
    """
    Find pseudos within top, return :class:`PseudoTable` object sorted by atomic number Z.
    """
    from monty.os.path import find_exts
    from pymatgen.io.abinitio.pseudos import PseudoTable, Pseudo
    exts=("psp8",)
    pseudos = []
    for p in find_exts(top, exts, exclude_dirs="_*"):
        try:
            pseudos.append(Pseudo.from_file(p))
        except Exception as exc:
            from warnings import warn
            warn("Exception in pseudo %s:\n%s" % (p.filepath, exc))
            
    return PseudoTable(pseudos).sort_by_z()


from collections import defaultdict, OrderedDict, Mapping
from monty.design_patterns import singleton
from pseudo_dojo.core.pseudos import DojoTable
from pseudo_dojo.pseudos import dojotable_absdir


class TableMetadata(object):
    """
     Metadata used to build one of the official tables.
    """
    def __init__(self, table_dir, djson_name):
        """
        Args:
            table_dir:
            djson_name:
        """
        self.table_dir = table_dir
        self.dojo_absdir = dojotable_absdir(table_dir)
        self.djson_name = djson_name
        self.djson_path = os.path.join(self.dojo_absdir, djson_name)

    #def __set__(self, obj, val):
    #    """Read-only data descriptor"""
    #    raise AttributeError("Dojo Tables are read-only!")

    #def __get__(self, obj, cls):
    #    print("in get with obj ", obj, "and cls", cls)
    #    x = cls
    #    #x = obj
    #    if not hasattr(x, "_dtables"):
    #        x._dtables = defaultdict(dict)

    #    if (self.table_dir in x._dtables and 
    #        self.djson_name in x._dtables[self.table_dir]):
    #        print("Returning cached table")
    #        return x._dtables[self.table_dir][self.djson_name]

    #    print("Not found, will compute and store the table")

    #    new_table = 1
    #    new_table = DojoTable.from_djson(self.djson_path)

    #    x._dtables[self.table_dir][self.djson_name] = new_table
    #    return new_table


@singleton
class OfficialTables(Mapping):

    def __init__(self):
        self._tables = tables = OrderedDict()

        # Naming scheme
        # [PP_TYPE]-[TABLE_NAME]-[XC_NAME]-[SO|""]-version
        #
        # where: 
        #
        #  PP_TYPE in ("PAW", "NC") defines the pseudos type.
        #  TABLE_NAME is the name of the table e.g. ONCVPSP, JTH ...
        #  XC_NAME gives the functional type e.g. PBE
        #  SO is present if the table contains the Spin-orbit term.
        #  maybe FR|SR for fully or scalar-relativistic?
        #  version is the version number.
        # 
        # Maybe descriptor?
        tables["ONC-DOJO-PBE"] = TableMetadata("ONCVPSP-PBE", "accuracy.djson")

        #tables["ONCV-GYGYv0.2-PBE"] = TableMetadata("ONCVPSP-PBE", "accuracy.djson")
        #tables["PAW-JTHv0.2-PBE"] = TableMetadata("ONCVPSP-PBE", "accuracy.djson")
        #tables["PAW-GBRVv0.2-PBE"] = TableMetadata("ONCVPSP-PBE", "accuracy.djson")

    # ABC protocol.
    def __iter__(self):
        return self._tables.__iter__()

    def __len__(self):
        return self._tables.__len__()

    def __getitem__(self, key):
        v = self._tables[key]
        if not isinstance(v, TableMetadata): return v 

        # Parse files, construct table and store it.
        new_table = DojoTable.from_djson(v.djson_path)
        new_table.dojo_name = key
        self._tables[key] = new_table
        return new_table

    #def nc_tables(self, **kwargs):
    #    return [table for for name, table in self.items() if table.dojo_info.isnc]

    #def paw_tables(self, **kwargs):
    #    return [table for for name, table in self.items() if table.dojo_info.ispaw]
