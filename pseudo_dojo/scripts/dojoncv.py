#!/usr/bin/env python
"""Script to generate/analyze/plot ONCVPSP pseudopotentials."""
import sys
import os
import collections
import argparse
import json
import shutil

from monty.termcolor import cprint
from abipy.flowtk.pseudos import Pseudo
from pseudo_dojo.core.dojoreport import DojoReport
from pseudo_dojo.ppcodes.ppgen import OncvGenerator
from pseudo_dojo.ppcodes.oncvpsp import OncvOutputParser, PseudoGenDataPlotter, oncv_make_open_notebook


def find_oncv_output(path):
    """
    Fix possible error in the specification of filename when we want a `.out` file.
    Return output path.
    """
    if path.endswith(".out"): return path
    root, _ = os.path.splitext(path)
    new_path = root + ".out"
    if not os.path.exists(new_path):
        raise ValueError("Cannot find neither %s nor %s" % (path, new_path))
    cprint("Maybe you meant %s" % new_path, "yellow")
    return new_path


def oncv_nbplot(options):
    """Generate jupyter notebook to plot data. Requires oncvpsp output file."""
    out_path = find_oncv_output(options.filename)
    return oncv_make_open_notebook(out_path)


def oncv_gnuplot(options):
    """Plot data with gnuplot."""
    out_path = find_oncv_output(options.filename)

    # Parse output file.
    onc_parser = OncvOutputParser(out_path)
    onc_parser.scan()
    if not onc_parser.run_completed:
        cprint("oncvpsp output is not complete. Exiting", "red")
        return 1

    onc_parser.gnuplot()
    return 0


def oncv_plot(options):
    """Plot data with matplotlib. Requires oncvpsp output file."""
    out_path = find_oncv_output(options.filename)

    # Parse output file.
    onc_parser = OncvOutputParser(out_path)
    onc_parser.scan()
    if not onc_parser.run_completed:
        cprint("oncvpsp output is not complete. Exiting", "red")
        return 1

    # Build the plotter
    plotter = onc_parser.make_plotter()

    # Plot data
    plotter.plot_radial_wfs()
    plotter.plot_atanlogder_econv()
    plotter.plot_projectors()
    plotter.plot_potentials()
    #plotter.plot_der_potentials()
    #for order in [1,2,3,4]:
    #    plotter.plot_der_densities(order=order)
    plotter.plot_densities()
    #plotter.plot_densities(timesr2=True)
    plotter.plot_den_formfact()
    return 0

    # Table of methods
    callables = collections.OrderedDict([
        ("wp", plotter.plot_waves_and_projs),
        ("dp", plotter.plot_dens_and_pots),
        ("lc", plotter.plot_atanlogder_econv),
        ("df", plotter.plot_den_formfact),
    ])

    # Call function depending on options.plot_mode
    if options.plot_mode == "slide":
        for func in callables.values():
            func()
    else:
        func = callables.get(options.plot_mode, None)
        if func is not None:
            func()
        else:
            plotter.plot_key(key=options.plot_mode)


def oncv_json(options):
    """
    Produce a string with the results in a JSON dictionary and exit
    Requires oncvpsp output file.
    """
    out_path = find_oncv_output(options.filename)
    onc_parser = OncvOutputParser(out_path)
    onc_parser.scan()
    if not onc_parser.run_completed:
        cprint("oncvpsp output is not complete. Exiting", "red")
        return 1

    # Generate json files with oncvpsp results.
    print(json.dumps(onc_parser.to_dict, indent=-1))
    return 0


def oncv_run(options):
    """
    Run oncvpsp, generate djrepo file, plot results. Requires input file.
    """
    # Select calc_type
    calc_type = dict(nor="non-relativistic",
                     sr="scalar-relativistic",
                     fr="fully-relativistic")[options.rel]

    # Build names of psp8 and djson files from input and relativistic mode.
    in_path = options.filename
    root, _ = os.path.splitext(in_path)

    # Enforce convention on output files.
    if options.rel == "nor":
        if not root.endswith("_nor"): root += "_nor"
    elif options.rel == "fr":
        if not root.endswith("_r"):
            root += "_r"
            cprint("FR calculation with input file without `_r` suffix. Will add `_r` to output files", "yellow")

    # Build names of output files.
    psp8_path = root + ".psp8"
    djrepo_path = root + ".djrepo"
    out_path = root + ".out"
    if os.path.exists(psp8_path):
        cprint("%s already exists and will be overwritten" % os.path.relpath(psp8_path), "yellow")
    if os.path.exists(djrepo_path):
        cprint("%s already exists and will be overwritten" % os.path.relpath(djrepo_path), "yellow")
    if os.path.exists(out_path):
        cprint("%s already exists and will be overwritten" % os.path.relpath(out_path), "yellow")

    # Build Generator and start generation.
    oncv_ppgen = OncvGenerator.from_file(in_path, calc_type, workdir=None)
    print(oncv_ppgen)
    print(oncv_ppgen.input_str)

    oncv_ppgen.start()
    retcode = oncv_ppgen.wait()

    if oncv_ppgen.status != oncv_ppgen.S_OK:
        cprint("oncvpsp returned %s. Exiting" % retcode, "red")
        return 1

    # Tranfer final output file.
    shutil.copy(oncv_ppgen.stdout_path, out_path)

    # Parse the output file
    onc_parser = OncvOutputParser(out_path)
    onc_parser.scan()
    if not onc_parser.run_completed:
        cprint("oncvpsp output is not complete. Exiting", "red")
        return 1

    # Extract psp8 files from the oncvpsp output and write it to file.
    s = onc_parser.get_psp8_str()
    with open(psp8_path, "wt") as fh:
        fh.write(s)

    # Write upf if available.
    upf_str = onc_parser.get_upf_str()
    if upf_str is not None:
        with open(psp8_path.replace(".psp8", ".upf"), "wt") as fh:
            fh.write(upf_str)

    pseudo = Pseudo.from_file(psp8_path)
    if pseudo is None:
        cprint("Cannot parse psp8 file: %s" % psp8_path, "red")
        return 1

    # Initialize and write djson file.
    report = DojoReport.empty_from_pseudo(pseudo, onc_parser.hints, devel=False)
    report.json_write()

    return 0


def main():

    def str_examples():
        return """\
Usage example:
    dojoncv.py run H.in         ==> Run oncvpsp input file (scalar relativistic mode).
    dojoncv.py plot H.out       ==> Use matplotlib to plot oncvpsp results for pseudo H.psp8.
    dojoncv.py gnuplot H.out    ==> Use gnuplot to plot oncvpsp results for pseudo H.psp8.
    dojoncv.py nbplot H.out     ==> Generate jupyter notebook to plot oncvpsp results.
    dojoncv.py json H.out       ==> Generate JSON file.
"""

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg: sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    # Parent parser implementing common options.
    copts_parser = argparse.ArgumentParser(add_help=False)
    copts_parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
                        help='Verbose, can be supplied multiple times to increase verbosity')

    copts_parser.add_argument('--loglevel', default="ERROR", type=str,
                        help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")

    copts_parser.add_argument('filename', default="", help="Path to the output file")

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=str_examples(), formatter_class=argparse.RawDescriptionHelpFormatter)

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    # Create the parsers for the sub-commands
    p_run = subparsers.add_parser('run', parents=[copts_parser], help=oncv_run.__doc__)
    p_run.add_argument("--rel", default="sr", help=("Relativistic treatment: `nor` for non-relativistic, "
                       "`sr` for scalar-relativistic, `fr` for fully-relativistic."))
    #p_run.add_argument("-d", "--devel", action="store_true", default=False,
    #                    help="put only two energies in the ecuts list for testing for developing the pseudo")

    # Create the parsers for the sub-commands
    p_plot = subparsers.add_parser('plot', parents=[copts_parser], help=oncv_plot.__doc__)
    p_plot.add_argument("-p", "--plot-mode", default="slide",
                        help=("Quantity to plot. Possible values: %s" %
                              str(["slide", "wp, dp, lc"] + PseudoGenDataPlotter.all_keys) + "\n"
                              "wp --> wavefunctions and projectors\n" +
                              "dp --> densities and potentials\n" +
                              "lc --> atan(logder) and convergence wrt ecut\n" +
                              "df --> density form factor"))

    p_nbplot = subparsers.add_parser('nbplot', parents=[copts_parser], help=oncv_nbplot.__doc__)

    p_gnuplot = subparsers.add_parser('gnuplot', parents=[copts_parser], help=oncv_gnuplot.__doc__)

    p_json = subparsers.add_parser('json', parents=[copts_parser], help=oncv_json.__doc__)

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

    # Dispatch
    return globals()["oncv_" + options.command](options)


if __name__ == "__main__":
    sys.exit(main())
