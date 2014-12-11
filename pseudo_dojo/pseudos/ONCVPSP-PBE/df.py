#!/usr/bin/env python
from __future__ import division, print_function, unicode_literals

import sys
import os
import pandas as pd

from tabulate import tabulate
from collections import OrderedDict, namedtuple
from pprint import pprint
from pymatgen.io.abinitio import Pseudo


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

    from pymatgen.io.abinitio.pseudos import PseudoTable
    pseudos = PseudoTable(pseudos)
    data, errors = pseudos.get_dojo_dataframe()
    print(data)

    if errors:
        print("ERRORS:")
        pprint(errors)

    accuracies = ["low", "normal", "high"]
    keys = ["dfact_meV", "v0", "b0_GPa", "b1", "ecut"]
    columns = ["symbol"] + [acc + "_" + k for k in keys for acc in accuracies]
    #print(columns)

    #rows, names, errors = [], [], []
    #for p in pseudos:
    #    report = p.read_dojo_report()
    #    df_entry = report.get("deltafactor", None)
    #    if df_entry is None:
    #        errors.append((p.name, "no deltafactor"))
    #        continue

    #    try:
    #        d = {"symbol": p.symbol}
    #        for acc in accuracies:
    #            d[acc + "_ecut"] = report["hints"][acc]["ecut"]

    #        for acc in accuracies:
    #            for k in keys:
    #                if k == "ecut": continue
    #                d[acc + "_" + k] = float(df_entry[acc][k])
    #        #print(d)
    #        names.append(p.name)
    #        rows.append(d)

    #    except Exception as exc:
    #        #raise
    #        print(p.name, "exc", str(exc))
    #        errors.append((p.name, str(exc)))

    ##print(rows)
    #data = pd.DataFrame(rows, index=names, columns=columns)
    #data = data[data["high_dfact_meV"] <= data["high_dfact_meV"].mean()]
    #data = data[data["high_dfact_meV"] <= 9]

    wrong = data[data["high_b1"] < 0]
    if not wrong.empty:
        print("WRONG".center(80, "*") + "\n", wrong)

    data = data[
        [acc + "_dfact_meV" for acc in accuracies]
        + [acc + "_ecut" for acc in accuracies]
    ]

    print("\nONCVPSP TABLE:\n") #.center(80, "="))
    columns = [acc + "_dfact_meV" for acc in accuracies] 
    columns += [acc + "_ecut" for acc in accuracies] 
    #print(data.to_string(columns=columns))
    tablefmt = "grid"
    print(tabulate(data[columns], headers="keys", tablefmt=tablefmt))

    print("\nSTATS:\n") #.center(80, "="))
    #print(data.describe())
    print(tabulate(data.describe(), headers="keys", tablefmt=tablefmt))

    bad = data[data["high_dfact_meV"] > data["high_dfact_meV"].mean()]
    print("\nPSEUDOS with high_dfact > mean:\n") # ".center(80, "*"))
    #print(bad)
    print(tabulate(bad, headers="keys", tablefmt=tablefmt))





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
