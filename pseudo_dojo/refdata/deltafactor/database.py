"""
This module provides a database with the deltafactor results,
and tools to compute the deltafactor of a pseudopotential.
Client code should use the official API df_database(xc) to access the database.

Example::
    db = df_database("PBE")
    wien2k = db.get_entry("Si")
    print(wien2k.v0, wien2k.b0, wien2k.b1)
"""
from __future__ import print_function, division, unicode_literals

import os
import os.path
import collections
import numpy as np

from monty.functools import lazy_property
from pymatgen.core.units import FloatWithUnit
from pymatgen.core.xcfunc import XcFunc
from pymatgen.util.plotting import add_fig_kwargs


class DeltaFactorEntry(collections.namedtuple("DeltaFactorEntry", "symbol v0 b0 b1 xc")):
    """
    Namedtuple storing the volume, the bulk-modulus and the pressure derivative b1.

    v0 is in A**3/natom, b0 is in eV /A**3, b1 is dimensionless.
    xc is the :class:`XcFunc` object associated to these data.
    """
    def __new__(cls, *args):
        """Extends the base class adding type conversion of arguments."""
        new_args = len(args) * [None]

        for i, arg in enumerate(args[:-1]):
            converter = float
            if i == 0: converter = str
            new_args[i] = converter(arg)

        v0 = FloatWithUnit(new_args[1], "ang^3")
        b0 = FloatWithUnit(new_args[2], "eV ang^-3")

        return super(cls, DeltaFactorEntry).__new__(cls,
                     symbol=new_args[0], v0=v0, b0=b0, b1=new_args[3], xc=args[-1])

    @property
    def b0_GPa(self):
        """b0 in GPa units."""
        return self.b0.to("GPa")

    @lazy_property
    def dfact_meV(self):
        """Deltafactor in meV"""
        return df_wien2k(self.symbol, self.v0, self.b0_GPa, self.b1, self.xc, b0_GPa=True)

    @property
    def dfactprime_meV(self):
        """Renormalized deltafactor proposed by Jollet."""
        return self.dfact_meV * (30 * 100) / (self.v0 * self.b0_GPa),


def read_data_from_filepath(filepath, xc):
    """
    Reads (v0, b0, b1) from file filename. b0 is in GPa.
    xc is the XcFunc associated to these data.
    Returns a dict of `DeltaFactorEntry` objects indexed by element symbol.
    """
    data = collections.OrderedDict()

    with open(filepath, "r") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("#") or not line:
                continue
            tokens = line.split()
            symbol = tokens[0]
            # Conversion GPa --> eV / A**3 and pass xc info to the Entry.
            tokens[2] = FloatWithUnit(tokens[2], "GPa").to("eV ang^-3")
            tokens.append(xc)
            data[symbol] = DeltaFactorEntry(*tokens)

    return data


def read_tables_from_file(filepath, xc):
    """
    Read deltafactor data from filepath. Returns pandas dataframe.
    """
    import pandas as pd
    columns = ["v0", "b0_GPa", "b1", "deltaf"]
    title = None
    tables = []

    with open(filepath, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line: continue
            #print(line)

            if line.startswith("#"):
                #if title is not None:
                #    table = pd.DataFrame(rows, index=symbols, columns=columns)
                #    table.title, table.info = title, info
                #    tables.append(table)

                title, info = line.split("--")
                symbols, rows = [], []
                continue

            tokens = line.split()
            # Conversion GPa --> eV / A**3
            symbol = tokens[0]
            symbols.append(symbol)
            b0_GPa = float(tokens[2])
            #b0 = FloatWithUnit(b0_GPa, "GPa").to("eV ang^-3")
            row = {"v0": float(tokens[1]), "b0_GPa": b0_GPa, "b1": float(tokens[3])}

            df = df_wien2k(symbol, row["v0"], b0_GPa, row["b1"], xc, b0_GPa=True)
            row.update({"deltaf": df})
            rows.append(row)

    table = pd.DataFrame(rows, index=symbols, columns=columns)
    table.title, table.info = title, info
    tables.append(table)

    for table in tables:
        print(table.title)
        print(table)

    return tables


class DeltaFactorDatabaseError(Exception):
    """Exceptions raised by the database."""


class DeltaFactorDatabase(object):
    """
    This object stores the deltafactor results for a given code and XC.
    It provides methods to access the data and plot the results.
    This object is not meant to instantiated explicitly.
    Client code should get the database via the public API:

        wien2k_db = df_database("PBE")
    """
    # Reference code.
    _REF_CODE = "WIEN2k"

    # mapping xc_name --> file_list
    _FILES4XC = {
        "PBE": ["WIEN2k-PBE.txt", "VASP-PBE.txt"],
        "PW": ["WIEN2k-PW.txt",],
    }

    # mapping xc_name --> location of CIF files.
    _CIFDIR4XC = {
        "PBE": "CIFs-PBE",
        "PW": "CIFs-PW",
    }

    # Supports other LDA flavors
    for xclda in ["PW_MOD",]:
        _FILES4XC[xclda] = _FILES4XC["PW"]
        _CIFDIR4XC[xclda] = _CIFDIR4XC["PW"]

    # Hack needed to support PBEsol.
    for xcgga in ["PBEsol",]:
        _FILES4XC[xcgga] = _FILES4XC["PBE"]
        _CIFDIR4XC[xcgga] = _CIFDIR4XC["PBE"]

    del xclda

    Error = DeltaFactorDatabaseError

    def __init__(self, xc):
        """"xc specifies the XC functionals used."""
        self.xc = XcFunc.asxc(xc)
        dirpath = os.path.abspath(os.path.dirname(__file__))
        self.dirpath = os.path.join(dirpath, "data")

        self._data = {}
        for bname in self._FILES4XC[self.xc]:
            file_path = os.path.join(self.dirpath, bname)
            if os.path.isfile(file_path) and file_path.endswith(".txt"):
                code, ext = os.path.splitext(bname)
                code = code.split("-")[0]
                self._data[code] = read_data_from_filepath(file_path, self.xc)

        self._cif_paths = {}

        #loc = "CIFs" + "-" + str(self.xc)
        loc = self._CIFDIR4XC[self.xc]
        cif_dirpath = os.path.join(self.dirpath, loc)
        for bname in os.listdir(cif_dirpath):
            if bname.endswith(".cif"):
                symbol, ext = os.path.splitext(bname)
                self._cif_paths[symbol] = os.path.join(cif_dirpath, bname)

    def has_symbol(self, symbol):
        """True if we have an entry for this symbol"""
        return symbol in self._data[self._REF_CODE]

    def get_cif_path(self, symbol):
        """Returns the path to the CIF file used for the given symbol."""
        return self._cif_paths[symbol]

    @property
    def symbols(self):
        """Chemical symbols available in the database."""
        return list(self._cif_paths.keys())

    @lazy_property
    def codes(self):
        """List of code names e.g. WIEN2k, VASP"""
        return list(sorted(self._data.keys()))

    def get_entry(self, symbol, code=None):
        """
        Return a :class:`DeltaFactorEntry` for the given chemical symbol.

        Args:
            symbol: Chemical symbol
            code: String identifying the code used to compute the entry.
                Default is self._REF_CODE (Wien2K)

        raise:
            self.Error if symbol is not in the database.
        """
        if code is None: code = self._REF_CODE
        try:
            return self._data[code][symbol]
        except KeyError:
            raise self.Error("No entry found for code %s, symbol %s" % (code, symbol))

    def spinat_spinmode_for_symbol(self, symbol):
        """
        Return the value of spinat and spinmode from the chemical symbol.
        """
        # TODO: Not sure spinat is ok if LDA.
        # Include spin polarization for O, Cr and Mn (antiferromagnetic)
        # and Fe, Co, and Ni (ferromagnetic).
        # antiferromagnetic Cr, O, ferrimagnetic Mn
        spinat = None
        spin_mode = "unpolarized"

        if symbol in ["Fe", "Co", "Ni"]:
            spin_mode = "polarized"
            if symbol == "Fe":
                spinat = 2 * [(0, 0, 2.3)]
            if symbol == "Co":
                spinat = 2 * [(0, 0, 1.2)]
            if symbol == "Ni":
                spinat = 4 * [(0, 0, 0.6)]

        if symbol in ["O", "Cr", "Mn"]:
            # Here we could have problems with include_so since we don't enforce "afm"
            spin_mode = "afm"
            if symbol == 'O':
                spinat = [(0, 0, 1.5), (0, 0, 1.5), (0, 0, -1.5), (0, 0, -1.5)]
            elif symbol == 'Cr':
                spinat = [(0, 0, 1.5), (0, 0, -1.5)]
            elif symbol == 'Mn':
                spinat = [(0, 0, 2.0), (0, 0, 1.9), (0, 0, -2.0), (0, 0, -1.9)]

        return spinat, spin_mode

    @add_fig_kwargs
    def plot_error_of_code(self, codename_or_data, values=("v0", "b0", "b1"), ref_code=None, **kwargs):
        """Plot the error of one code vs ref_code."""
        import matplotlib.pyplot as plt

        if ref_code is None: ref_code = self._REF_CODE
        ref_data = self._data[ref_code]

        data = codename_or_data
        if isinstance(codename_or_data, str):
            data = self._data[codename_or_data]

        entries = list(ref_data.values())

        # Build grid of plots.
        fig, ax_list = plt.subplots(nrows=len(values), ncols=1, sharex=True, squeeze=False)
        ax_list = ax_list.ravel()

        for aname, ax in zip(values, ax_list):
            # Sort entries according to the value of the attribute aname.
            entries.sort(key=lambda t: getattr(t, aname))

            ord_symbols = [r.symbol for r in entries]
            xticks = []

            for i, osym in enumerate(ord_symbols):
                ref_value = getattr(ref_data[osym], aname)
                value = getattr(data[osym], aname)
                err = 100 * (value - ref_value) / ref_value
                ax.plot(i, err, "ro")
                xticks.append(i)

            # Set xticks and labels.
            ax.plot(xticks, np.zeros(len(xticks)), "b-")
            ax.set_xticks(xticks)
            ax.set_xticklabels(ord_symbols, rotation="vertical")

            #ax.grid(True)
            ax.set_ylabel("Relative error %s" % aname)
            #ax.set_xlabel("Ecut [Ha]")

        return fig


def check_cif_wien2k_consistency(xc):
    """
    This function compares the v0 of the structures read from the CIF files with
    with those computed from the text files containing the deltafactor results.
    """
    from pymatgen.core.structure import Structure
    import pandas as pd
    pd.set_option("display.max_rows", None)

    rows, index = [], []

    db = df_database(xc)
    for symbol in db.symbols:
        if symbol == "S": continue
        wien2k = db.get_entry(symbol)
        #vasp = db.get_entry(symbol, "VASP")
        vasp = db.get_entry(symbol, "VASP-relaxed")
        structure = Structure.from_file(db.get_cif_path(symbol))
        v0 = structure.volume / len(structure)
        wien2k_v0abs_diff = wien2k.v0 - v0
        wien2k_v0rel_diff = 100 * (wien2k.v0 - v0) / v0
        #if abs_(v0diff):
        print("Wien2k Symbol: %s, abs: %.2f, rel %.2f" % (symbol, wien2k_v0abs_diff, wien2k_v0rel_diff))

        vasp_v0abs_diff = vasp.v0 - v0
        vasp_v0rel_diff = 100 * (vasp.v0 - v0) / v0
        print("Vasp Symbol: %s, abs: %.2f, rel %.2f" % (symbol, vasp_v0abs_diff, vasp_v0rel_diff))
        #if symbol == "O":
        #    print(v0, wien2k.v0)
        index.append(symbol)
        rows.append(dict(wien2k_v0abs_diff=wien2k_v0abs_diff,
                         wien2k_v0rel_diff=wien2k_v0rel_diff,
                         vasp_v0abs_diff=vasp_v0abs_diff,
                         vasp_v0rel_diff=vasp_v0rel_diff,
                    ))

    frame = pd.DataFrame(data=rows, index=index)
    print(frame)
    print(frame.describe())

    import matplotlib.pyplot as plt
    frame.hist()
    plt.show()


##########################################################################################
# Official API to access the database.
##########################################################################################

# Mapping XC --> Database.
__DELTAF_DATABASE_XC = None


def df_database(xc):
    """
    Returns the deltafactor database with the reference results associated to this XC functional
    xc is the exchange-correlation functional e.g. PBE, PW
    """
    global __DELTAF_DATABASE_XC
    if __DELTAF_DATABASE_XC is None:
        __DELTAF_DATABASE_XC = {}

    # Create xc database and cache it.
    if xc not in __DELTAF_DATABASE_XC:
        db = DeltaFactorDatabase(xc=xc)
        __DELTAF_DATABASE_XC[db.xc] = db

    return __DELTAF_DATABASE_XC[xc]


def df_wien2k(symbol, v0f, b0f, b1f, xc, b0_GPa=False, v=3, useasymm=False):
    """
    Compute the deltafactor wrt to WIEN2k
    xc specifies the XC functional to be used.
    """
    wien2k = df_database(xc).get_entry(symbol)
    b0 = wien2k.b0
    if b0_GPa: b0 = wien2k.b0_GPa
    #print(v0f, b0f, b1f, wien2k.v0, b0, wien2k.b1)

    return df_compute(float(wien2k.v0), float(b0), wien2k.b1, float(v0f), b0f, b1f, b0_GPa=b0_GPa,
                      v=v, useasymm=useasymm)


def df_compute(v0w, b0w, b1w, v0f, b0f, b1f, b0_GPa=False, v=3, useasymm=False):
    """
    Compute the deltafactor. Based on the code of the offical calcDelta.py script.

    Args:
        v0w, b0w, b1w:  Volume, bulk-modulus and pressure derivative of b0w (reference values).
        v0f, b0f, b1f: Volume, bulk-modulus and pressure derivative of b0f (computed values).
        v: version of delta factor, current version is 3, 1 for symmetrical old version

    .. note::

        v0 is A**3/natom, by default b0 is in eV/A**3, GPa units are used if b0_GPa is True.
    """
    #print(v0w, b0w, b1w, v0f, b0f, b1f)

    if v == 1:
        # delta factor form version 1
        if b0_GPa:
            # Conversion GPa --> eV/A**3
            b0w = FloatWithUnit(b0w, "GPa").to("eV ang^-3")
            b0f = FloatWithUnit(b0f, "GPa").to("eV ang^-3")

        Vi = 0.94 * v0w
        Vf = 1.06 * v0w

        a3f = 9. * v0f**3. * b0f / 16. * (b1f - 4.)
        a2f = 9. * v0f**(7./3.) * b0f / 16. * (14. - 3. * b1f)
        a1f = 9. * v0f**(5./3.) * b0f / 16. * (3. * b1f - 16.)
        a0f = 9. * v0f * b0f / 16. * (6. - b1f)

        a3w = 9. * v0w**3. * b0w / 16. * (b1w - 4.)
        a2w = 9. * v0w**(7./3.) * b0w / 16. * (14. - 3. * b1w)
        a1w = 9. * v0w**(5./3.) * b0w / 16. * (3. * b1w - 16.)
        a0w = 9. * v0w * b0w / 16. * (6. - b1w)

        x = [0, 0, 0, 0, 0, 0, 0]

        x[0] = (a0f - a0w)**2
        x[1] = 6. * (a1f - a1w) * (a0f - a0w)
        x[2] = -3. * (2. * (a2f - a2w) * (a0f - a0w) + (a1f - a1w)**2.)
        x[3] = -2. * (a3f - a3w) * (a0f - a0w) - 2. * (a2f - a2w) * (a1f - a1w)
        x[4] = -3./5. * (2. * (a3f - a3w) * (a1f - a1w) + (a2f - a2w)**2.)
        x[5] = -6./7. * (a3f - a3w) * (a2f - a2w)
        x[6] = -1./3. * (a3f - a3w)**2.

        Fi = np.zeros_like(Vi)
        Ff = np.zeros_like(Vf)

        for n in range(7):
            Fi = Fi + x[n] * Vi**(-(2.*n-3.)/3.)
            Ff = Ff + x[n] * Vf**(-(2.*n-3.)/3.)

        return 1000. * np.sqrt((Ff - Fi) / (0.12 * v0w))

    elif v == 3:
        # version 3

        if b0_GPa:
            # Conversion GPa --> eV/A**3
            b0w = FloatWithUnit(b0w, "GPa").to("eV ang^-3")
            b0f = FloatWithUnit(b0f, "GPa").to("eV ang^-3")

        if useasymm:
            Vi = 0.94 * v0w
            Vf = 1.06 * v0w
        else:
            Vi = 0.94 * (v0w + v0f) / 2.
            Vf = 1.06 * (v0w + v0f) / 2.

        a3f = 9. * v0f**3. * b0f / 16. * (b1f - 4.)
        a2f = 9. * v0f**(7./3.) * b0f / 16. * (14. - 3. * b1f)
        a1f = 9. * v0f**(5./3.) * b0f / 16. * (3. * b1f - 16.)
        a0f = 9. * v0f * b0f / 16. * (6. - b1f)

        a3w = 9. * v0w**3. * b0w / 16. * (b1w - 4.)
        a2w = 9. * v0w**(7./3.) * b0w / 16. * (14. - 3. * b1w)
        a1w = 9. * v0w**(5./3.) * b0w / 16. * (3. * b1w - 16.)
        a0w = 9. * v0w * b0w / 16. * (6. - b1w)

        x = [0, 0, 0, 0, 0, 0, 0]

        x[0] = (a0f - a0w)**2
        x[1] = 6. * (a1f - a1w) * (a0f - a0w)
        x[2] = -3. * (2. * (a2f - a2w) * (a0f - a0w) + (a1f - a1w)**2.)
        x[3] = -2. * (a3f - a3w) * (a0f - a0w) - 2. * (a2f - a2w) * (a1f - a1w)
        x[4] = -3./5. * (2. * (a3f - a3w) * (a1f - a1w) + (a2f - a2w)**2.)
        x[5] = -6./7. * (a3f - a3w) * (a2f - a2w)
        x[6] = -1./3. * (a3f - a3w)**2.

        y = [0, 0, 0, 0, 0, 0, 0]

        y[0] = (a0f + a0w)**2 / 4.
        y[1] = 3. * (a1f + a1w) * (a0f + a0w) / 2.
        y[2] = -3. * (2. * (a2f + a2w) * (a0f + a0w) + (a1f + a1w)**2.) / 4.
        y[3] = -(a3f + a3w) * (a0f + a0w) / 2. - (a2f + a2w) * (a1f + a1w) / 2.
        y[4] = -3./20. * (2. * (a3f + a3w) * (a1f + a1w) + (a2f + a2w)**2.)
        y[5] = -3./14. * (a3f + a3w) * (a2f + a2w)
        y[6] = -1./12. * (a3f + a3w)**2.

        Fi = np.zeros_like(Vi)
        Ff = np.zeros_like(Vf)

        Gi = np.zeros_like(Vi)
        Gf = np.zeros_like(Vf)

        for n in range(7):
            Fi += x[n] * Vi**(-(2.*n-3.)/3.)
            Ff += x[n] * Vf**(-(2.*n-3.)/3.)

            Gi += y[n] * Vi**(-(2.*n-3.)/3.)
            Gf += y[n] * Vf**(-(2.*n-3.)/3.)

        Delta = 1000. * np.sqrt((Ff - Fi) / (Vf - Vi))

        #Deltarel = 100. * np.sqrt((Ff - Fi) / (Gf - Gi))
        #if useasymm:
        #    Delta1 = 1000. * np.sqrt((Ff - Fi) / (Vf - Vi)) \
        #             / v0w / b0w * vref * bref
        #else:
        #    Delta1 = 1000. * np.sqrt((Ff - Fi) / (Vf - Vi)) \
        #             / (v0w + v0f) / (b0w + b0f) * 4. * vref * bref

        return Delta

    else:
        raise ValueError("Wrong version %s" % v)


if __name__ == "__main__":
    #import sys
    #read_tables_from_file(sys.argv[1], xc)
    check_cif_wien2k_consistency()
