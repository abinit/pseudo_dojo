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
