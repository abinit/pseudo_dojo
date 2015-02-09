from __future__ import division, print_function, unicode_literals
#from .dojo import Dojo

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
