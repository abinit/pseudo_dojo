# coding: utf-8
"""Public API to access the pseudopotential tables."""
from __future__ import division, print_function, unicode_literals
# TODO
import os

from pymatgen.core.periodic_table import PeriodicTable
from pymatgen.io.abinitio.pseudos import PseudoTable #, Pseudo



class DojoTables(object):
    """This object gathers the different tables in a single namespace."""


_PSEUDOS_EXTS = set(("psp8",))


class DojoPseudoTable(PseudoTable):
    """
    A pseudopotential table provided by pseudo_dojo.
    We subclass PseudoTable so that we can easily add extra methods, if needed.
    """


def get_nc_table(nc_type, xc_type, table_type=""):
    top = "-".join([nc_type, xc_type])
    if table_type: top = "-".join([top, table_type])
    top = os.path.join(os.path.dirname(__file__), top)
    if not os.path.isdir(top):
        raise ValueError("Directory %s does not exist")

    # Gather all pseudos starting from the current working directory 
    table = PeriodicTable()
    all_symbols = set(element.symbol for element in table.all_elements)
    dirs = [os.path.join(top, d) for d in os.listdir(top) if d in all_symbols]

    # directory: find all pseudos with the psp8 extensions ignore directories starting with _
    #top = paths[0]
    #paths, ext = [], "psp8"
    #for dirpath, dirnames, filenames in os.walk(top):
    #    if os.path.basename(dirpath).startswith("_"): continue
    #    dirpath = os.path.abspath(dirpath)
    #    for filename in filenames:
    #        if any(filename.endswith(ext) for ext in exts):
    #            paths.append(os.path.join(dirpath, filename))

    paths = []
    for dir in dirs:
        paths.extend(os.path.join(dir, f) for f in os.listdir(dir) if any(f.endswith(ext) for ext in _PSEUDOS_EXTS))

    #print(paths)
    return DojoPseudoTable(paths).sort_by_z()


if __name__ == "__main__":
    print(get_nc_table("ONCVPSP", "PBE"))
