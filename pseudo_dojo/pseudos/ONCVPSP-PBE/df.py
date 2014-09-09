#!/usr/bin/env python
from __future__ import print_function, division

import sys
import os
import pandas as pd

from collections import OrderedDict, namedtuple
from pprint import pprint
from pymatgen.io.abinitio import Pseudo


#class DojoEntry(object):
#    def __init_(self, path):
#        self.path = os.path.abspath(path)
#        self.name  = os.path.basename(path)
#        self.table = os.path.basename(os.path.dirname(self.path))
#        self.pseudo = Pseudo.from_file(self.path) 


def main():
    try:
        top = sys.argv[1]
    except IndexError:
        top = "."

    pseudos = []
    for dirpath, dirnames, filenames in os.walk(top):
        # Exclude pseudos in _inputs
        if os.path.basename(dirpath) == "_inputs": continue

        #print(filenames)
        pseudos.extend([Pseudo.from_file(os.path.join(dirpath, f)) for f in filenames 
            if f.endswith(".psp8")])
            #if f.endswith(".psp8") and "-" not in f])
            #if f.endswith(".psp8") and "-" in f])
    #print(pseudos)

    #accuracies = ["low", "normal", "high"]
    accuracies = ["normal", "high"]
    keys = ["dfact_meV", "v0", "b0_GPa", "b1", "ecut"]
    columns = ["symbol"] + [acc + "_" + k for k in keys for acc in accuracies]
    #print(columns)

    rows, names, errors = [], [], []
    for p in pseudos:
        report = p.read_dojo_report()
        df_entry = report.get("deltafactor", None)
        if df_entry is None:
            errors.append((p.name, "no deltafactor"))
            continue

        try:
            d = {"symbol": p.symbol}
            for acc in accuracies:
                d[acc + "_ecut"] = report["hints"][acc]["ecut"]

            for acc in accuracies:
                for k in keys:
                    if k == "ecut": continue
                    d[acc + "_" + k] = float(df_entry[acc][k])
            #print(d)
            names.append(p.name)
            rows.append(d)

        except Exception as exc:
            #raise
            print(p.name, "exc", str(exc))
            errors.append((p.name, str(exc)))

    #print(rows)
    data = pd.DataFrame(rows, index=names, columns=columns)
    #data = data[data["high_dfact_meV"] <= data["high_dfact_meV"].mean()]
    #data = data[data["high_dfact_meV"] <= 9]

    wrong = data[data["high_b1"] < 0]
    if not wrong.empty:
        print("WRONG".center(80, "*") + "\n", wrong)

    data = data[
        [acc + "_dfact_meV" for acc in accuracies]
        + [acc + "_ecut" for acc in accuracies]
    ]

    print("ALL DATA".center(80, "="))
    columns = [acc + "_dfact_meV" for acc in accuracies] 
    columns += [acc + "_ecut" for acc in accuracies] 

    print(data.to_string(columns=columns))
    print(data.describe())

    bad = data[data["high_dfact_meV"] > data["high_dfact_meV"].mean()]
    print("BAD".center(80, "*"))
    print(bad)

    if errors:
        print("ERRORS:")
        pprint(errors)

    #import matplotlib.pyplot as plt 
    #bad.plot(kind="barh")
    #bad.plot(kind="kde")
    #bad.plot(kind="density")
    #bad["high_dfact_meV"].plot(kind="bar")
    ###bad[[acc + "_dfact_meV" for acc in ["normal", "high"]]].plot(kind="bar")
    #bad.plot(kind="bar", columns=["high_dfact_meV"])
    #bad["delta_normal"].hist(bins=200)
    #bad["delta_high"].hist(bins=200)
    #data["high_dfact_meV"].hist(bins=200)
    #plt.show()


if __name__ == "__main__":
    sys.exit(main())
