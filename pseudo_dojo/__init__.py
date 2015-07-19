from __future__ import division, print_function, unicode_literals


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


class OfficialTable(object):
    """
    A data descriptor that sets and returns values
    normally and prints a message logging their access.
    """
    def __init__(self, initval=None, name='var'):
        self.val = initval
        self.name = name

    def __get__(self, obj, objtype):
        print("obj ", obj, "objtype", objtype)
        print('Retrieving', self.name)
        return self.val

    def __set__(self, obj, val):
        """Eead-only data descriptor"""
        #print('Updating', self.name)
        #self.val = val
        raise AttributeError("Dojo Tables are read-only!")


class Tables(object):
    """
    This object gathers the official tables provided by PseudoDojo in a single namespace.
    """
    GGA = OfficialTable(initval="hello")


