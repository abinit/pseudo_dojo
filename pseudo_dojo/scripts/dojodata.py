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
from monty.os.path import find_exts
from pymatgen.io.abinitio.pseudos import Pseudo
from pseudo_dojo.core.pseudos import DojoTable


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
            #report.plot_etotal_vs_ecut(title=pseudo.basename)
            if options.eos: 
                report.plot_deltafactor_eos(title=pseudo.basename)
            #fig = report.plot_deltafactor_convergence(title=pseudo.basename, what="dfactprime_meV", show=True)
            fig = report.plot_deltafactor_convergence(title=pseudo.basename, show=True)
            #report.plot_deltafactor_convergence(title=pseudo.basename)

        #from bokeh import mpl
        #from bokeh.plotting import show
        #show(mpl.to_bokeh(fig, name="test"))
        #return

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

        # phonon
        if any(k in options.what_plot for k in ("all", "phonon")):
            if report.has_trial("phonon"):
                report.plot_phonon_convergence(title=pseudo.basename)

        #if any(k in options.what_plot for k in ("all", "phwoa")):
        #    if report.has_trial("phwoa"):
        #        report.plot_phonon_convergence(title=pseudo.basename+"W/O ASR", woasr=True)



def dojo_compare(options):
    """Plot and compare DOJO results for multiple pseudos."""
    pseudos = options.pseudos
    for z in pseudos.zlist: 
        pseudos_z = pseudos[z]
        if len(pseudos_z) > 1:
            pseudos_z.dojo_compare(what=options.what_plot)
        else:
            print("Find only one pseudo for Z=%s" % z)


def dojo_trials(options):
    """Visualize the results of the different tests."""
    pseudos = options.pseudos

    # Build pandas DataFrame
    data, errors = pseudos.get_dojo_dataframe()
    #print(data)

    if errors:
        print("ERRORS:")
        pprint(errors)
    
    #import matplotlib.pyplot as plt
    #data.plot_trials(savefig=options.savefig)
    #data.plot_hist(savefig=options.savefig)
    #data.sns_plot(savefig=options.savefig)
    print(data["high_dfact_meV"])
    import matplotlib.pyplot as plt
    ax = data["high_dfact_meV"].hist()
    ax.set_xlabel("Deltafactor [meV]")
    plt.show()


def dojo_table(options):
    """Build and show a pandas table."""
    pseudos = options.pseudos

    data, errors = pseudos.get_dojo_dataframe()
    #data.tabulate()
    #return

    if False:
        """Select best entries"""
        #best = {}
        #symbols = set(data["symbol"])
        #for sym in symbols:
        grouped = data.groupby("symbol")
        #print(grouped.groups)

        rows, names = [], []
        for name, group in grouped:
            #print(name, group["high_dfact_meV"])
            best = group.sort("high_dfact_meV").iloc[0]
            names.append(name)
            #print(best.keys())
            #print(best.name)
            l = {k: getattr(best, k) for k in ("name", "Z", 'high_b0_GPa', 'high_b1', 'high_v0', 'high_dfact_meV', 'high_ecut')} 
            rows.append(l)

        import pandas
        best_frame = pandas.DataFrame(rows, index=names)
        best_frame = best_frame.sort("Z")
        print(tabulate(best_frame, headers="keys"))
        print(tabulate(best_frame.describe(),  headers="keys"))

        import matplotlib.pyplot as plt
        best_frame["high_dfact_meV"].hist(bins=100)
        plt.show()

        return

    if errors:
        print("ERRORS:")
        pprint(errors)

    accuracies = ["low", "normal", "high"]
    accuracies = ["low"]
    keys = ["dfact_meV", "v0", "b0_GPa", "b1", "ecut"]
    columns = ["symbol"] + [acc + "_" + k for k in keys for acc in accuracies]
    print(columns)

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

    accuracies = ["low", "high"]
    data = data[
        [acc + "_dfact_meV" for acc in accuracies]
      + [acc + "_ecut" for acc in accuracies]
#      + [acc + "_gbrv_fcc_a0_rel_err" for acc in accuracies]
#      + [acc + "_gbrv_bcc_a0_rel_err" for acc in accuracies]
    ]

    print("\nONCVPSP TABLE:\n") #.center(80, "="))
    columns = [acc + "_dfact_meV" for acc in accuracies] 
    columns += [acc + "_ecut" for acc in accuracies] 
#    columns += [acc + "_gbrv_fcc_a0_rel_err" for acc in accuracies] 
#    columns += [acc + "_gbrv_bcc_a0_rel_err" for acc in accuracies] 

    tablefmt = "grid"
    floatfmt=".2f"
    print(tabulate(data[columns], headers="keys", tablefmt=tablefmt, floatfmt=floatfmt))
    #print(data.to_string(columns=columns))

    print("\nSTATS:\n") #.center(80, "="))
    #print(data.describe())
    print(tabulate(data.describe(), headers="keys", tablefmt=tablefmt, floatfmt=floatfmt))

    bad = data[data["high_dfact_meV"] > data["high_dfact_meV"].mean()]
    print("\nPSEUDOS with high_dfact > mean:\n") # ".center(80, "*"))
    print(tabulate(bad, headers="keys", tablefmt=tablefmt, floatfmt=floatfmt))

    #gbrv_fcc_bad = data[data["high_gbrv_fcc_a0_rerr"] > (data["high_gbrv_fcc_a0_rerr"].abs()).mean()]
    #print("\nPSEUDOS with high_dfact > mean:\n") # ".center(80, "*"))
    #print(tabulate(bad, headers="keys", tablefmt=tablefmt, floatfmt=floatfmt))


def dojo_dist(options):
    """
    Plot the distribution of the deltafactor and of the relative error for the GBRV fcc/bcc tests.
    """
    fig = options.pseudos.plot_dfgbrv_dist()


def dojo_check(options):

    for p in options.pseudos:
        try:
            report = p.dojo_report
        except Exception as exc:
            print("Invalid dojo_report in:", p.basename)
            print("Exception: ", exc)
            continue

        # Comment this to fix the md5 checksum in the pseudos
        #p.check_and_fix_dojo_md5()

        #if "ppgen_hints" not in report: # and "deltafactor" not in report:
        #    print(p.basename, "old version without ppgen_hints")
        #    continue

        try:
            error = report.check()
            if error:
                print("[%s] Validation problem" % p.basename)
                print(error)
                print()

        except Exception as exc:
            print("Error: ", p.basename + str(exc))


def dojo_validate(options):
    for pseudo in options.pseudos:
        try:
            report = pseudo.dojo_report
            #if report.is_validated:
            #    print("Pseudo %s is already validated!" % pseudo.basename)
            #report.plot_deltafactor_convergence(title=pseudo.basename, what="dfactprime_meV")

            hints = report.compute_hints()
            print("hints for %s computed from deltafactor prime: %s" % (pseudo.basename, hints))

            #ans = prompt("Do you accept the hints? [Y]")
            #ans = False
            #if ans:
            #    print("got true")
            #    #report.validate(hints)
            #    #pseudo.write_dojo_report(report)
            #else:
            #    print("The dojoreport contains ecuts :\n%s" % report.ecuts)
            #    #new_ecuts = prompt("Enter new ecuts to compute (comma-separated values or empty string to abort)")
            #    #if new_ecuts:
            #    #    print("Exit requested by user")
            #    #    return 

            #    # Be careful with the format here! it should be %.1f
            #    #new_ecuts = np.array(new_ecuts.split(","))
            #    #report.add_ecuts(new_ecuts)
            #    #pseudo.write_dojo_report(report)

        except Exception as exc:
            print(pseudo.basename, "raised: ", str(exc))


def main():
    def str_examples():
        examples = """\
Usage example:\n
    dojodata plot H.psp8                ==> Plot dojo data for pseudo H.psp8
    dojodata trials H.psp8 -r 1
    dojodata compare H.psp8 H-low.psp8  ==> Plot and compare dojo data for pseudos H.psp8 and H-low.psp8
    dojodata table .                    ==> Build table (find all psp8 files withing current directory)
"""
        return examples

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg: sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    # Parent parser for commands that need to know on which subset of pseudos we have to operate.
    def parse_rows(s):
        if not s: return []
        tokens = s.split(",")
        return list(map(int, tokens)) if tokens else []

    def parse_symbols(s):
        if not s: return []
        return s.split(",")

    pseudos_selector_parser = argparse.ArgumentParser(add_help=False)
    pseudos_selector_parser.add_argument('pseudos', nargs="+", help="Pseudopotential file or directory containing pseudos")
    pseudos_selector_parser.add_argument('-s', "--symbols", type=parse_symbols, help=("List of chemical symbols to include or exclude."
        "Example --symbols=He,Li to include He and Li, --symbols=-He to exclude He"))

    # Options for pseudo selection.
    group = pseudos_selector_parser.add_mutually_exclusive_group()
    group.add_argument("-r", '--rows', default="", type=parse_rows, help="Select these rows of the periodic table.")
    group.add_argument("-f", '--family', type=str, default="", help="Select this family of the periodic table.")

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=str_examples(), formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--loglevel', default="ERROR", type=str,
                        help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")

    parser.add_argument('--seaborn', action="store_true", help="Use seaborn settings")

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    plot_options_parser = argparse.ArgumentParser(add_help=False)
    plot_options_parser.add_argument("-w", "--what-plot", type=str, default="all", 
                                      help="Quantity to plot e.g df for deltafactor, gbrv for GBRV tests")
    plot_options_parser.add_argument("-e", "--eos", action="store_true", help="Plot EOS curve")

    # Subparser for plot command.
    p_plot = subparsers.add_parser('plot', parents=[pseudos_selector_parser, plot_options_parser], help="Plot DOJO_REPORT data.")

    # Subparser for compare.
    p_compare = subparsers.add_parser('compare', parents=[pseudos_selector_parser, plot_options_parser], help="Compare pseudos")

    # Subparser for table command.
    p_table = subparsers.add_parser('table', parents=[pseudos_selector_parser], help="Build pandas table.")

    # Subparser for dist command.
    p_dist = subparsers.add_parser('dist', parents=[pseudos_selector_parser], 
                                   help="Plot distribution of deltafactor and GBRV relative errors.")

    # Subparser for trials command.
    p_trials = subparsers.add_parser('trials', parents=[pseudos_selector_parser], help="Plot DOJO trials.")
    p_trials.add_argument("--savefig", type=str, default="", help="Save plot to savefig file")

    # Subparser for check command.
    p_check = subparsers.add_parser('check', parents=[pseudos_selector_parser], help="Check pseudos")

    # Subparser for validate command.
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

    def get_pseudos(options):
        """
        Find pseudos in paths, return :class:`DojoTable` object sorted by atomic number Z.
        Accepts filepaths or directory.
        """
        exts=("psp8",)

        paths = options.pseudos
        if len(paths) == 1 and os.path.isdir(paths[0]):
            top = paths[0]
            paths = find_exts(top, exts, exclude_dirs="_*")
            #table = DojoTable.from_dir(paths[0])

        pseudos = []
        for p in paths:
            try:
                pseudos.append(Pseudo.from_file(p))
            except Exception as exc:
                warn("Error in %s:\n%s. This pseudo will be ignored" % (p, exc))

        table = DojoTable(pseudos)

        # Here we select a subset of pseudos according to family or rows
        if options.rows:
            table = table.select_rows(options.rows)
        elif options.family:
            table = table.select_families(options.family)

        if options.symbols:
            table = table.select_symbols(options.symbols)

        return table.sort_by_z()

    # Build DojoTable from the paths specified by the user.
    options.pseudos = get_pseudos(options)

    if options.seaborn:
        import seaborn as sns
        #sns.set(style='ticks', palette='Set2')
        sns.set(style="dark", palette="Set2")
        #And to remove "chartjunk", do:
        #sns.despine()
        #plt.tight_layout()
        #sns.despine(offset=10, trim=True)

    # Dispatch
    globals()["dojo_" + options.command](options)
    return 0


if __name__ == "__main__":
    sys.exit(main())
