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


from collections import defaultdict
from pseudo_dojo.core.pseudos import DojoTable
from pseudo_dojo.pseudos import dojotable_absdir


class OfficialTable(object):
    """
    A data descriptor that sets and returns values
    normally and prints a message logging their access.
    """
    def __init__(self, table_dir, djson_name):
        self.table_dir = table_dir
        self.dojo_absdir = dojotable_absdir(table_dir)
        self.djson_name = djson_name
        self.djson_path = os.path.join(self.dojo_absdir, djson_name)

    def __set__(self, obj, val):
        """Read-only data descriptor"""
        raise AttributeError("Dojo Tables are read-only!")

    def __get__(self, obj, cls):
        print("in get with obj ", obj, "and cls", cls)
        x = cls
        #x = obj
        if not hasattr(x, "_dtables"):
            x._dtables = defaultdict(dict)

        if (self.table_dir in x._dtables and 
            self.djson_name in x._dtables[self.table_dir]):
            print("Returning cached table")
            return x._dtables[self.table_dir][self.djson_name]

        print("Not found, will compute and store the table")

        new_table = 1
        new_table = DojoTable.from_djson(self.djson_path)

        x._dtables[self.table_dir][self.djson_name] = new_table
        return new_table


class Tables(object):
    """
    This object gathers the official tables provided by PseudoDojo in a single namespace.
    """
    GGA = OfficialTable("ONCVPSP-PBE", "accuracy.djson")

    #NC_ODRH_GGA_SR_V0_2 = OfficialTable("ONCVPSP-PBE", "accuracy.djson")
    #PAW_JTH_PBE_V0_2

    #@classmethod
    #def objects(cls):
    #    dtables = cls._dtables

    #    all_tables = []
    #    for v in dtables.values():
    #        all_tables.extend(list(v.values()))

    #    return all_tables

    #@classmethod
    #def get_all(cls)

