from __future__ import division, print_function

import sys
import os
import os.path
import collections
import numpy as np

##########################################################################################

class DeltaFactorEntry(collections.namedtuple("DeltaFactorEntry", "symbol v0 b0 bp")):
    "Namedtuple storing the volume, the bulk-modulus and Bp for given symbol (a.u.)"

    def __new__(cls, *args):
        "Extends the base class adding type conversion of arguments."
        new_args = len(args) * [None]

        for (i, arg) in enumerate(args):
            converter = float
            if i == 0: converter = str
            new_args[i] = converter(arg)

        return super(cls, DeltaFactorEntry).__new__(cls, *new_args)

def read_data_from_filename(filename):
    """
    Reads (v0, b0, bp) from file filename
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
            data[symbol] = DeltaFactorEntry(*tokens)
    return data


def singleton(cls):
    """This decorator can be used to create a singleton out of a class."""
    instances = {}
    def getinstance():
        if cls not in instances:
            instances[cls] = cls()
        return instances[cls]
    return getinstance

##########################################################################################

@singleton
class DeltaFactorDataset(object):
    """This object stores the deltafactor results."""
    # Reference code.
    _REF_CODE = "WIEN2k"

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
                #print(codename)
                d[code] = read_data_from_filename(file_path)

        self.cif_paths = d = {}

        cif_dirpath = os.path.join(self.dirpath, "CIFs")
        for entry in os.listdir(cif_dirpath):
            if entry.endswith(".cif"):
                symbol, ext = os.path.splitext(entry)
                d[symbol] = os.path.join(cif_dirpath, entry)

    def get_cif_path(self, symbol):
        return self.cif_paths[symbol]

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
        return self._data[code][symbol]

    def plot_error_of_code(self, codename_or_data, values=("v0", "b0", "bp"), ref_code=None, **kwargs):
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
        #print(entries)

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

def compute_deltaf(v0w, b0w, b1w, v0f, b0f, b1f):
    # GPa --> SI
    b0w *= 10.**9. / 1.602176565e-19 / 10.**30.
    b0f *= 10.**9. / 1.602176565e-19 / 10.**30.

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

    Delta1 = 1000. * np.sqrt((Ff - Fi) / (0.12 * v0w))
    return Delta1
