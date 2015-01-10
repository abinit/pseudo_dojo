#!/usr/bin/env python
from __future__ import division, print_function, unicode_literals

import sys
import os
import argparse
import pandas as pd

from warnings import warn
from collections import OrderedDict, namedtuple
from pprint import pprint
from tabulate import tabulate
from pymatgen.io.abinitio.pseudos import PseudoTable 


def dojo_plot(options):
    top = options.path
    if not os.path.isdir(top):
        paths = [top]
    else:
        # directory: find all pseudos with the psp8 extensions.
        # ignore directories starting with _
        paths, ext = [], ".psp8"
        for dirpath, dirnames, filenames in os.walk(top):
            if os.path.basename(dirpath).startswith("_"): continue
            dirpath = os.path.abspath(dirpath)
            for filename in filenames:
                if filename.endswith(ext): 
                    paths.append(os.path.join(dirpath, filename))

    pseudos = PseudoTable(paths).sort_by_z()

    for pseudo in pseudos:
        if not pseudo.has_dojo_report:
            warn("pseudo %s does not contain the DOJO_REPORT section" % pseudo.filepath)
            continue

        report = pseudo.dojo_report
        # FIXME add symbol
        report.symbol = pseudo.symbol
        print(pseudo)
        print(report)

        print("trials passed: ", report.trials)
        report.has_hints
        #report.has_exceptions
        report.print_table()

        if False and report.has_trial("deltafactor"):
            report.plot_deltafactor_eos(title=pseudo.basename)
            report.plot_deltafactor_convergence(title=pseudo.basename)
            #report.plot_etotal_vs_ecut()

        for struct_type in ("fcc", "bcc"):
            if report.has_trial("gbrv_" + struct_type):
                report.plot_gbrv_eos(struct_type=struct_type, title=pseudo.basename)


def dojo_table(options):
    top = options.path
    
    if not os.path.isdir(top):
        pseudos = [top]
    else:
        pseudos = []
        for dirpath, dirnames, filenames in os.walk(top):
            # Exclude pseudos in _inputs
            if os.path.basename(dirpath) == "_inputs": continue
            pseudos.extend([os.path.join(dirpath, f) for f in filenames 
                if f.endswith(".psp8")])
                #if f.endswith(".psp8") and "-" not in f])
                #if f.endswith(".psp8") and "-" in f])

    #print(pseudos)
    pseudos = PseudoTable(pseudos).sort_by_z()
    data, errors = pseudos.get_dojo_dataframe()
    print(data)

    if errors:
        print("ERRORS:")
        pprint(errors)

    accuracies = ["low", "normal", "high"]
    keys = ["dfact_meV", "v0", "b0_GPa", "b1", "ecut"]
    columns = ["symbol"] + [acc + "_" + k for k in keys for acc in accuracies]
    #print(columns)

    #data = data[data["high_dfact_meV"] <= data["high_dfact_meV"].mean()]
    #data = data[data["high_dfact_meV"] <= 9]

    def calc_rerrors(data):
        # Relative error
        data["low_dfact_abserr"] = data["low_dfact_meV"] - data["high_dfact_meV"]
        data["normal_dfact_abserr"] =  data["normal_dfact_meV"] - data["high_dfact_meV"]
        data["low_dfact_rerr"] = 100 * (data["low_dfact_meV"] - data["high_dfact_meV"]) / data["high_dfact_meV"]
        data["normal_dfact_rerr"] = 100 * (data["normal_dfact_meV"] - data["high_dfact_meV"]) / data["high_dfact_meV"]

        for k in ["v0", "b0_GPa", "b1"]:
            data["low_" + k + "_abserr"] = data["low_" + k] - data["high_" + k]
            data["normal_" + k + "_abserr"] = data["normal_" + k] - data["high_" + k]
            data["low_" + k + "_rerr"] = 100 * (data["low_" + k] - data["high_" + k]) / data["high_" + k]
            data["normal_" + k + "_rerr"] = 100 * (data["normal_" + k] - data["high_" + k]) / data["high_" + k]

        return data

    import matplotlib.pyplot as plt
    #import seaborn as sns
    #data.plot(x="symbol", y="high_dfact_meV", use_index=True)
    #data = calc_rerrors(data)
    #g = sns.PairGrid(data, x_vars="Z", y_vars=[
    #    "low_ecut",
    #    "low_dfact_meV",
    #    #"normal_ecut",
    #    #"low_dfact_meV",
    #    #"high_dfact_meV", 
    #    "low_v0_rerr", "low_b0_GPa_rerr", "low_b1_rerr",
    #    ]
    #) #, hue="smoker")
    #g.map(plt.scatter)
    #g.add_legend()
    #plt.show()

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


def main():
    def str_examples():
        examples = """\
Usage example:\n
"""
        return examples

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg: sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    # Parent parser for commands that need to know on which subset of pseudos we have to operate.
    pseudos_selector_parser = argparse.ArgumentParser(add_help=False)
    group = pseudos_selector_parser.add_mutually_exclusive_group()
    #group.add_argument("-n", '--nids', default=None, type=parse_nids, help=(
    #    "Node identifier(s) used to select the task. Integer or comma-separated list of integers. Use `status` command to get the node ids.\n"
    #    "Examples: --nids=12 --nids=12,13,16 --nids=10:12 to select 10 and 11, --nids=2:5:2 to select 2,4"  
    #    ))

    group.add_argument('path', nargs="?", help="Pseudopotential file or directory containing pseudos")

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=str_examples(), formatter_class=argparse.RawDescriptionHelpFormatter)

    #parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
    #                    help='verbose, can be supplied multiple times to increase verbosity')

    parser.add_argument('--loglevel', default="ERROR", type=str,
                        help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    # Subparser for plot command.
    p_plot = subparsers.add_parser('plot', parents=[pseudos_selector_parser], help="Plot data.")

    # Subparser for table command.
    p_plot = subparsers.add_parser('table', parents=[pseudos_selector_parser], help="Build pandas table.")

    # Parse command line.
    try:
        options = parser.parse_args()
    except Exception as exc: 
        show_examples_and_exit(error_code=1)

    # loglevel is bound to the string value obtained from the command line argument. 
    # Convert to upper case to allow the user to specify --loglevel=DEBUG or --loglevel=debug
    import logging
    numeric_level = getattr(logging, options.loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % options.loglevel)
    logging.basicConfig(level=numeric_level)

    if options.command == "plot":
        dojo_plot(options)

    elif options.command == "table":
        dojo_table(options)

    else:
        raise ValueError("Don't know how to handle command %s" % options.command)

    return 0


if __name__ == "__main__":
    sys.exit(main())
