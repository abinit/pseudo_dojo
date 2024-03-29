#!/usr/bin/env python
from __future__ import division, print_function, unicode_literals

import sys
import os

from tabulate import tabulate
from pprint import pprint
from abipy.flowtk import Pseudo


def main():
    top = "."
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
    keys = ["dfact_meV", "v0", "b0_GPa", "b1", "ecut", "fcc_a0_rel_err", "bcc_a0_rel_err"]
    columns = ["symbol"] + [acc + "_" + k for k in keys for acc in accuracies]
    #print(columns)

    ##print(rows)
    #data = pd.DataFrame(rows, index=names, columns=columns)
    #data = data[data["high_dfact_meV"] <= data["high_dfact_meV"].mean()]
    #data = data[data["high_dfact_meV"] <= 9]

    def calc_rerrors(data):
        # Relative error
        data["low_dfact_abserr"] = data["low_dfact_meV"] - data["high_dfact_meV"]
        data["normal_dfact_abserr"] = data["normal_dfact_meV"] - data["high_dfact_meV"]
        data["low_dfact_rerr"] = 100 * (data["low_dfact_meV"] - data["high_dfact_meV"]) / data["high_dfact_meV"]
        data["normal_dfact_rerr"] = 100 * (data["normal_dfact_meV"] - data["high_dfact_meV"]) / data["high_dfact_meV"]

        for k in ["v0", "b0_GPa", "b1"]:
            data["low_" + k + "_abserr"] = data["low_" + k] - data["high_" + k]
            data["normal_" + k + "_abserr"] = data["normal_" + k] - data["high_" + k]
            data["low_" + k + "_rerr"] = 100 * (data["low_" + k] - data["high_" + k]) / data["high_" + k]
            data["normal_" + k + "_rerr"] = 100 * (data["normal_" + k] - data["high_" + k]) / data["high_" + k]

        return data

    #import seaborn as sns
    import matplotlib.pyplot as plt
    #data = calc_rerrors(data)
    #g = sns.PairGrid(data, x_vars="Z", y_vars=[
    #    "low_ecut",
    #    "low_dfact_meV",
    #    #"normal_ecut",
    #    #"low_dfact_meV",
    #    #"high_dfact_meV",
    #    #"low_v0_rerr", "low_b0_GPa_rerr", "low_b1_rerr",
    #    ]
    #) #, hue="smoker")
    #g.map(plt.scatter)
    #g.add_legend()

    #data["high_dfact_meV"].hist(bins=200)
    #data["high_fcc_a0_rel_err"].hist(bins=200)
    #data["high_bcc_a0_rel_err"].hist(bins=200)

    keys = [
        "dfact_meV",
        #"dfactprime_meV",
        "bcc_a0_rel_err", "fcc_a0_rel_err",
        #"ecut",
    ]

    fig, ax_list = plt.subplots(nrows=len(keys), ncols=1, sharex=False, squeeze=False)
    ax_list = ax_list.ravel()

    #kind = "density"
    kind = "scatter"
    kind = "line"
    zmin, zmax = 1, 110
    for ax, key in zip(ax_list, keys):
        for acc in accuracies:
            c = data[acc + "_" + key][data.Z <= zmax]
            c = c[data.Z >= zmin]
            c.plot(kind=kind, ax=ax, style="o-", legend=True)
            #c.hist(ax=ax, bins=200)

    plt.show()

    wrong = data[data["high_b1"] < 0]
    if not wrong.empty:
        print("WRONG".center(80, "*") + "\n", wrong)

    data = data[
        [acc + "_dfact_meV" for acc in accuracies]
        + [acc + "_ecut" for acc in accuracies]
        #+ [acc + "_fcc_a0_rel_err" for acc in accuracies]
        #+ [acc + "_bcc_a0_rel_err" for acc in accuracies]
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
