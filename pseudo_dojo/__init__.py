from __future__ import division, print_function, unicode_literals

import os

from collections import OrderedDict
from collections import abc
from monty.dev import deprecated
from monty.functools import lazy_property
from monty.design_patterns import singleton
#from pymatgen.core.periodic_table import Element
from pseudo_dojo.core.pseudos import OfficialDojoTable
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
class OfficialTables(abc.Mapping):
    """
    Official tables provided by the PseudoDojo project.

    Naming scheme:

        [PP_TYPE]-[XC_NAME]-[TABLE_NAME_WITH_VERSION]-[DJSON_NAME]

    where:

        #. PP_TYPE in ("PAW", "NC") defines the pseudos type.
        #. XC_NAME defines the exchange-correlation functional e.g. PBE
        #. TABLE_NAME_WITH_VERSION is the name of the table e.g. ONCVPSP, JTH ...
        #. DJSON_NAME is the name of the djson file located in the pseudodojo directory.

    In the case of NC pseudos, the SOC term is treated via djson e.g. "accuracy-soc.djson".
    maybe FR|SR for fully or scalar-relativistic?
    """
    def __init__(self):
        self._tables = tables = OrderedDict()

        #tables["ONCVPSP-PBE-PDv0.2-accuracy"] = TableMetadata("ONCVPSP-PBE-PDv0.2", "test_accuracy.djson")
        tables["ONCVPSP-PBE-PDv0.2-accuracy"] = TableMetadata("ONCVPSP-PBE-PDv0.3", "standard.djson")
        #tables["ONCV-GYGYv0.2-PBE"] = TableMetadata("ONC-PBE-GYGYv0.2", "htc.djson")
        #tables["PAW-JTHv0.2-PBE"] = TableMetadata("PAW-PBE-JTHv0.2", "standard.djson")
        #tables["PAW-PBE-GBRVv0.2"] = TableMetadata("PAW-PBE-GBRVv0.2", "efficiency.djson")
        # TODO: Add check on md5 of djson files.

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

        new_table = OfficialDojoTable.from_djson_file(v.djson_path)
        new_table.dojo_name = key
        self._tables[key] = new_table

        return new_table

    def select_tables(self, pp_type=None, xc=None):
        """
        Low-level method used to select tables.

        Args:
            pp_type: NC for norm-conserving pseudos or PAW. If None no selection is done.
            xc: If xc is not None, only the pseudos with this xc functional are selected.
                else no selection is performed.
        Return:
            List of tables.
        """
        tables = self.values()
        if pp_type is not None:
            tables = [t for t in tables if t.dojo_info.pp_type == pp_type]
        if xc is not None:
            tables = [t for t in tables if t.dojo_info.xc == xc]

        return tables

    @lazy_property
    def available_xcs(self):
        """List with the XC functionals available."""
        return sorted(set(t.xc for t in self.values()))

    def all_nctables(self, xc=None):
        """
        Return all norm-conserving tables available.
        If xc is not None, only the pseudos with this xc functional are selected.
        """
        tables = [t for t in self.values() if t.dojo_info.isnc]
        if xc is not None:
            tables = [t for t in tables if t.xc == xc]
        return tables

    def all_pawtables(self, xc=None):
        """
        Return all PAW tables available.
        If xc is not None, only the pseudos with this xc functional are selected.
        """
        tables = [table for table in self.values() if table.dojo_info.ispaw]
        if xc is not None:
            tables = [t for t in tables if t.xc == xc]
        return tables

    def pseudo_from_symbol_md5(self, symbol, md5, pp_type, xc):
        """
        Find pseudo from chemical symbol, md5 checksum, pp_type and xc.

        Raises:
            ValueError if pseudo is not found.
        """
        # This is the __eq__ implemented for Pseudo
        #return (self.md5 == other.md5 and
        #        self.__class__ == other.__class__ and
        #        self.Z == other.Z and
        #        self.Z_val == other.Z_val and
        #        self.l_max == other.l_max )

        # Here it's very important the we don't have duplicated md5 in the pseudo dojo tables.
        # Actually there's a check in test/test_dojotables.py
        for p in self.select_pseudos(symbol, pp_type=pp_type, xc=xc):
            if p.md5 == md5: return p

        raise ValueError(
            "Cannot find pseudo associated to the following parameters:\n"
            "symbol=%s, md5=%s, pp_type=%s, xc=%s" % (symbol, md5, pp_type, xc))

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

    def select_nc_pseudos(self, symbol, xc=None):
        """
        Return list of NC pseudos with the given chemical `symbol`.
        If xc is not None, only the pseudos with this xc functional are selected.
        """
        return self.select_pseudos(symbol, pp_type="NC", xc=xc)

    def select_paw_pseudos(self, symbol, xc=None):
        """
        Return list of PAW pseudos with the given chemical `symbol`.
        If xc is not None, only the pseudos with this xc are selected.
        """
        return self.select_pseudos(symbol, pp_type="PAW", xc=xc)

    #def plot_numpseudos(self, **kwargs):
    #    """Plot periodic table with the number of pseudos available per element."""
    #    import matplotlib.pyplot as plt
    #    from ptplotter.plotter import ElementDataPlotter

    #    symbols, data, max_num = [], {}, 0
    #    for el in Element:
    #        pseudos = self.select_pseudos(el.symbol, pp_type=None, xc=None)
    #        symbols.append(el.symbol)
    #        max_num = max(max_num, len(pseudos))
    #        data[el.symbol] = {"num_pseudos": len(pseudos)}

    #    def count_pseudos(elt):
    #        """Pseudos available"""
    #        return elt.get('num_pseudos', 0)

    #    epd = ElementDataPlotter(elements=symbols, data=data)
    #    #epd.ptable(functions=[count_pseudos], cmaps=["jet"], )
    #    #epd.cmaps[0].set_clim([0, max_num])
    #    #epd.redraw_ptable()

    #    plt.show()


def logo1():
    return """\
````````````````````````````````````````````````````````````````````````````````````````````````````
```````````````-shyo.```````````````````````````````````````````````````````````````````````````````
``````````-/oo`+hhhh-:s+/.````````````-:::::::-.````````````````````````````````````````-:::::.`````
```````.+ss/-.``-//-``.:+ys/.`````````yhhooooyhy+``......```.......`..````..``......````/oooo+-`````
`````.+yo-````````````````:sy:````````yhh`````shh.+ysooss+`+ysoooo+.yy````yy.:yyoooso-`/sso+oyo-````
````.sy-````````````````````/h+```````yhh::::/yhs.oyo+///:`+hs////..hh.```hh.:hs```-hy-hh.```/hs````
````sy-```````-------.```````/h/``````yhhooooo+:.`//::/ohy.+ho::::.`yh-..-hh.:hs...:hy.yh-``.+ho````
``-shysssssss+yhhhhhhosssssssshhs`````oss`````````/soooos+`/ssoooo+`-oooooo:`:ssoooo+.`-osooos+.````
``-hhhhhhhhhhshhhhhhhshhhhhhhhhhy`````````````````````````````````````````````````````````...```````
``-hhhhhhhhhhshhhhhhhshhhhhhhhhhy`````://////:-.```.//////`````````.////+/``````````````````````````
```-yy-------.+++++++:-------:hs-`````yhh++++shy+``.+ooo+/`````--``.+ooo+/``````````````````````````
````:ho`````````````````````.yy.``````yhh`````ohh--sy+//oy+````yh.-sy+//oyo`````````````````````````
`````-ys-``````````````````/yo.```````yhh`````+hh:sh/````yh-```yh-oh/````yh:````````````````````````
``````./ys/.````````````-/ys:`````````yhh:::/+hho`+ho.``:hh-..-hh.+ho.``-hh.````````````````````````
`````````:+sso/::-::/+oys+-```````````+oooooo+/-```:+oooo/.:ooo+-``:+oooo+.`````````````````````````
````````````.-//++++/:-````````````````````````````````````````````````````.```````````````.````````
"""


def logo2():
    return """\
MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMd/--+NMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMNdyo+No----dh+oymNMMMMMMMMMMMdhhhhhhhdNMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMdhhhhhNMMMMM
MMMMMMMms//ydNMMdyymMMmhs/+yNMMMMMMMMM/--++++/--sMMNNNNNNMMNNNNNNNNMNNMMMMNNMMNNNNNNMMMMyoo+oodMMMMM
MMMMMNs:omMMMMMMMMMMMMMMMNh/:hMMMMMMMM/--MMMMM/--mo:+oo+/sMs-/oooooN::NMMN:-Nh-:o++/+dMy//oo+/+mMMMM
MMMMN/:dMMMMMMMMMMMMMMMMMMMNy-sMMMMMMM/--hhhhs--/No:+osyyhMo-/yyyyNN--NMMN--Nh-/MMMd-:m--NMMMy-/MMMM
MMMN/-mNNNNNNMmmmmmmmmNNNNNNNy-yNMMMMM/--+++++ohNMyyhhyo-:No-+hhhhmN--dNNm--Nh-/NNNh-/m--dNMNs-+MMMM
MMm/--///////s-------+////////--+MMMMMo//MMMMMMMMMy++ooo/sNy/+oooooMd++oo++dMd//oo++smMdo+oo++sNMMMM
MMd----------+-------+-----------MMMMMMNMMMMMMMMMMMMNNNNMMMMNNNNNNNMMMNNNNMMMMNNNNMMMMMMMNNNNMMMMMMM
MMd----------+-------+-----------MMMMMhsssssyhdNMMMmsssssyMMMMMMMMMNsssssyMMMMMMMMMMMMMMMMMMMMMMMMMM
MMNd-:dmmmmmmmssssssshmmmmmmmh-+mMMMMM/--ssso/-:sMMmsooooyMMMMMmdMMNsooooyMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMh-oNMMMMMMMMMMMMMMMMMMMMm/:NMMMMMM/--MMMMNo--dd/:syy+:oMMMM:-md/:syy+:oNMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMd:/dMMMMMMMMMMMMMMMMMNy:+NMMMMMMM/--MMMMMs--h+-yMMMN:-dMMM--m+-yMMMM/-dMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMNy:/yNMMMMMMMMMMMMms:/hMMMMMMMMM/--hhhyo--+Ns-+mMNd--mmNd--Ns-+mMNd--mMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMNho//oyhhddhys+:/smMMMMMMMMMMMs+++++ooymMMMho+++osmhooosmMMho+++osmMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMNdyssossyhmNMMMMMMMMMMMMMMNMMMMNMMMMMMMMNMMMMMMMMMMMMMMMMMNMMMMNMMMMMMMMMMMMMMMNMMMMMMMM
"""
