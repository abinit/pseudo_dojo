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
from pymatgen.io.abinitio.pseudos import PseudoTable, Pseudo


def dojo_plot(options):
    """Plot DOJO results for a single pseudo."""
    pseudos = options.pseudos
    for pseudo in pseudos:
        #print(pseudos)
        if not pseudo.has_dojo_report:
            warn("pseudo %s does not contain the DOJO_REPORT section" % pseudo.filepath)
            continue

        report = pseudo.dojo_report
        #print(pseudo)
        #print(report)

        #print("trials passed: ", report.trials)
        #report.has_hints
        #report.has_exceptions
        #report.print_table()

        # Deltafactor
        if report.has_trial("deltafactor") and any(k in options.what_plot for k in ("all", "df")):
            report.plot_etotal_vs_ecut(title=pseudo.basename)
            if options.eos:
                report.plot_deltafactor_eos(title=pseudo.basename)
            report.plot_deltafactor_convergence(title=pseudo.basename)

        # GBRV
        if any(k in options.what_plot for k in ("all", "gbrv")):
            count = 0
            for struct_type in ("fcc", "bcc"):
                trial = "gbrv_" + struct_type
                if report.has_trial(trial):
                    count += 1
                    if options.eos:
                        report.plot_gbrv_eos(struct_type=struct_type, title=pseudo.basename)
            if count:
                report.plot_gbrv_convergence(title=pseudo.basename)


def dojo_compare_plots(options):
    """Plot and compare DOJO results for multiple pseudos."""
    pseudos = options.pseudos
    import matplotlib.pyplot as plt

    # Compare ecut convergence and Deltafactor
    if all(p.dojo_report.has_trial("deltafactor") for p in pseudos) and \
           any(k in options.what_plot for k in ("all", "df")):

        fig, ax_list = plt.subplots(nrows=len(pseudos), ncols=1, sharex=True, squeeze=True)
        for ax, pseudo in zip(ax_list, pseudos):
            pseudo.dojo_report.plot_etotal_vs_ecut(ax=ax, show=False, label=pseudo.basename)
        plt.show()

        fig, ax_grid = plt.subplots(nrows=5, ncols=len(pseudos), sharex=True, sharey="row", squeeze=False)
        for ax_list, pseudo in zip(ax_grid.T, pseudos):
            pseudo.dojo_report.plot_deltafactor_convergence(ax_list=ax_list, show=False)

        fig.suptitle(" vs ".join(p.basename for p in pseudos))
        plt.show()

    # Compare GBRV results
    if all(p.dojo_report.has_trial("gbrv_bcc") for p in pseudos) and \
       any(k in options.what_plot for k in ("all", "gbrv")):

        fig, ax_grid = plt.subplots(nrows=2, ncols=len(pseudos), sharex=True, sharey="row", squeeze=False)
        for ax_list, pseudo in zip(ax_grid.T, pseudos):
            pseudo.dojo_report.plot_gbrv_convergence(ax_list=ax_list, show=False)

        fig.suptitle(" vs ".join(p.basename for p in pseudos))
        plt.show()


def dojo_table(options):
    """Build table."""
    pseudos = options.pseudos

    # Build pandas DataFrame
    data, errors = pseudos.get_dojo_dataframe()
    print(data)

    import matplotlib.pyplot as plt
    #data.show_hist()
    #data.show_trials()
    #data.sns_plot()
    #plt.show()
    #return

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


def dojo_validate(options):
    errors = []
    for p in options.pseudos:
        report = p.dojo_report
        if "ppgen_hints" not in report: # and "deltafactor" not in report:
            print(p.basename, "old version")
            continue

        try:
            d = p.dojo_report.validate()
            if d:
                print("Validation problem -->", p.basename, d)
        except Exception as exc:
            print("Error: ", p.basename + str(exc))

    if errors:
        print(errors)


def main():
    def str_examples():
        examples = """\
Usage example:\n
    dojodata plot H.psp8             ==> Plot dojo data for pseudo H.psp8
    dojodata plot H.psp8 H-low.psp8  ==> Plot and compare dojo data for pseudos H.psp8 and H-low.psp8
    dojodata table .                 ==> Build table
"""
        return examples

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg: sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    # Parent parser for commands that need to know on which subset of pseudos we have to operate.
    pseudos_selector_parser = argparse.ArgumentParser(add_help=False)
    pseudos_selector_parser.add_argument('pseudos', nargs="+", help="Pseudopotential file or directory containing pseudos")

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=str_examples(), formatter_class=argparse.RawDescriptionHelpFormatter)

    #parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
    #                    help='verbose, can be supplied multiple times to increase verbosity')

    parser.add_argument('--loglevel', default="ERROR", type=str,
                        help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    # Subparser for plot command.
    p_plot = subparsers.add_parser('plot', parents=[pseudos_selector_parser], help="Plot DOJO_REPORT data.")
    p_plot.add_argument("-w", "--what-plot", type=str, default="all", help="Quantity to plot e.g df for deltafactor, gbrv for GBRV tests")
    p_plot.add_argument("-e", "--eos", type=bool, default=False, help="Plot EOS curve")

    # Subparser for table command.
    p_table = subparsers.add_parser('table', parents=[pseudos_selector_parser], help="Build pandas table.")

    # Subparser for table validate.
    p_validate = subparsers.add_parser('validate', parents=[pseudos_selector_parser], help="Validate pseudos")

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

    def read_pseudos(paths, exts):
        """
        Find pseudos in paths, return PseudoTable object sorted by Z.
        Accepts filepaths or directory.
        """
        if len(paths) == 1 and os.path.isdir(paths[0]):
            # directory: find all pseudos with the psp8 extensions.
            # ignore directories starting with _
            top = paths[0]
            paths, ext = [], "psp8"
            for dirpath, dirnames, filenames in os.walk(top):
                if os.path.basename(dirpath).startswith("_"): continue
                dirpath = os.path.abspath(dirpath)
                for filename in filenames:
                    if any(filename.endswith(ext) for ext in exts):
                        paths.append(os.path.join(dirpath, filename))

        pseudos = []
        for p in paths:
            try:
                pseudos.append(Pseudo.from_file(p))
            except:
                print("Error in %s" % p)
                                                                       
        return PseudoTable(pseudos).sort_by_z()

    # Build PseudoTable from the paths specified by the user.
    options.pseudos = read_pseudos(options.pseudos, exts=("psp8",))

    if options.command == "plot":
        if len(options.pseudos) > 1:
            dojo_compare_plots(options)
        else:
            dojo_plot(options)

    elif options.command == "table":
        dojo_table(options)

    elif options.command == "validate":
        dojo_validate(options)

    else:
        raise ValueError("Don't know how to handle command %s" % options.command)

    return 0


if __name__ == "__main__":
    sys.exit(main())
