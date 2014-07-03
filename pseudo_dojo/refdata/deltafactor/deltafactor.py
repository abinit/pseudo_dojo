"""
This module provides a databased for accessing the deltafactor results,
and tools to compute the deltafactor of a pseudopotential.
Client code should use the official API df_database() to access the database.

Example::
    db = df_database()
    wien2k = db.get_entry("Si")
    print(wien2k.v0, wien2k.b0, wien2k.b1)
"""
from __future__ import division, print_function

import sys
import os
import os.path
import collections
import numpy as np

from pymatgen.core.units import FloatWithUnit


class DeltaFactorEntry(collections.namedtuple("DeltaFactorEntry", "symbol v0 b0 b1")):
    """
    Namedtuple storing the volume, the bulk-modulus and the pressure derivative b1.

    v0 is in A**3/natom, b0 is in eV /A**3, b1 is dimensionless.
    """
    def __new__(cls, *args):
        """Extends the base class adding type conversion of arguments."""
        new_args = len(args) * [None]

        for (i, arg) in enumerate(args):
            converter = float
            if i == 0: converter = str
            new_args[i] = converter(arg)

        v0 = FloatWithUnit(new_args[1], "ang^3")
        b0 = FloatWithUnit(new_args[2], "eV ang^-3")

        new = super(cls, DeltaFactorEntry).__new__(cls, symbol=new_args[0], v0=v0, b0=b0, b1=new_args[3])
        return new

    @property
    def b0_GPa(self):
        """b0 in GPa units."""
        return self.b0.to("GPa") 

def read_data_from_filename(filename):
    """
    Reads (v0, b0, b1) from file filename. b0 is in GPa.
    Returns a dict of `DeltaFactorEntry` objects indexed by element symbol.
    """
    data = collections.OrderedDict()

    with open(filename, "r") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("#") or not line: 
                continue
            tokens = line.split()
            symbol = tokens[0] 
            # Conversion GPa --> eV / A**3
            tokens[2] = FloatWithUnit(tokens[2], "GPa").to("eV ang^-3") 
            data[symbol] = DeltaFactorEntry(*tokens)

    return data



class DeltaFactorDatabaseError(Exception):
    """Exceptions raised by the database."""


class DeltaFactorDatabase(object):
    """
    This object stores the deltafactor results and 
    provides methods to access the data.
    """
    # Reference code.
    _REF_CODE = "WIEN2k"

    Error = DeltaFactorDatabaseError

    def __init__(self):
        self.dirpath = os.path.abspath(os.path.dirname(__file__))
        self.dirpath = os.path.join(self.dirpath, "data")

        self._data = d = {}
        for entry in os.listdir(self.dirpath):
           file_path = os.path.join(self.dirpath, entry)
           if os.path.isfile(file_path) and file_path.endswith(".txt"):
                code, ext = os.path.splitext(entry)
                if code == "README": 
                    continue
                d[code] = read_data_from_filename(file_path)

        self._cif_paths = d = {}

        cif_dirpath = os.path.join(self.dirpath, "CIFs")
        for entry in os.listdir(cif_dirpath):
            if entry.endswith(".cif"):
                symbol, ext = os.path.splitext(entry)
                d[symbol] = os.path.join(cif_dirpath, entry)

    def has_symbol(self, symbol):
        """True if we have an entry for this symbol"""
        return symbol in self._data[self._REF_CODE]

    def get_cif_path(self, symbol):
        """Returns the path to the CIF file used for the given symbol."""
        return self._cif_paths[symbol]

    def get_entry(self, symbol, code=None):
        """
        Return the `DeltaFactorEntry` for the given chemical symbol.

        Args:
            symbol:
                Chemical symbol
            code:
                String identifying the code used to compute the entry. 
                Default is self._REF_CODE (Wien2K)
        """
        if code is None:
            code = self._REF_CODE 
        try:
            return self._data[code][symbol]
        except KeyError:
            raise self.Error("No entry found for code %s, symbol %s" % (code, symbol))

    def plot_error_of_code(self, codename_or_data, values=("v0", "b0", "b1"), ref_code=None, **kwargs):
        import matplotlib.pyplot as plt

        # Extract keyword arguments.
        show = kwargs.pop("show", True)
        savefig = kwargs.pop("savefig",None)
        title = kwargs.pop("title", None)

        if ref_code is None:
            ref_code = self._REF_CODE 

        ref_data = self._data[ref_code]

        data = codename_or_data
        if isinstance(codename_or_data, str):
            data = self._data[codename_or_data]

        entries = ref_data.values()

        # Build grid of plots.
        fig, ax_list = plt.subplots(nrows=len(values), ncols=1, sharex=True, squeeze=False)
        ax_list = ax_list.ravel()

        if title:
            fig.suptitle(title)

        for (aname, ax) in zip(values, ax_list):
            # Sort entries according to the value of the attribute aname.
            entries.sort(key = lambda t: getattr(t, aname))
            #print(entries)

            ord_symbols = [r.symbol for r in entries]
            xticks = []

            for (i, osym) in enumerate(ord_symbols):
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

        if show:
            plt.show()

        if savefig is not None:
            fig.savefig(savefig)
                         
        return fig

##########################################################################################
# Official API to access the database.
##########################################################################################

_DELTAF_DATABASE = DeltaFactorDatabase()


def df_database():
    """Returns the deltafactor database with the reference results."""
    return _DELTAF_DATABASE


def df_compute(v0w, b0w, b1w, v0f, b0f, b1f, b0_GPa=False, v=3, useasymm=False):
    """
    Compute the deltafactor. Based on the code of the offical calcDelta.py script.

    Args:
        v0w, b0w, b1w: 
            Volume, bulk-modulus and pressure derivative of b0w (reference values).
        v0f, b0f, b1f:
            Volume, bulk-modulus and pressure derivative of b0w (computed values).

    .. note:

        v0 is A**3/natom, by default b0 is in eV/A**3, GPa units are used if b0_GPa is True.
    """

    if v == 1:
        # delta factor form verion 1
        if b0_GPa:
            # Conversion GPa --> eV/A**3
            b0w = FloatWithUnit(b0w, "GPa").to("eV Ang^-3")
            b0f = FloatWithUnit(b0f, "GPa").to("eV Ang^-3")

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