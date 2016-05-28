#!/usr/bin/env python
"""Script to generate/analyze/plot ONCVPSP pseudopotentials."""
from __future__ import division, print_function, unicode_literals

import sys
import collections
import argparse
import json

from monty.termcolor import cprint
from monty.os.path import which
from pseudo_dojo.core.pseudos import dojopseudo_from_file
from pseudo_dojo.ppcodes.oncvpsp import OncvOutputParser, PseudoGenDataPlotter


def oncv_plot(options):
    """Plot data. Requires oncvpsp output file."""
    out_path = options.filename

    onc_parser = OncvOutputParser(out_path)
    onc_parser.scan()
    if not onc_parser.run_completed:
        cprint("oncvpsp output is not complete. Exiting", "red")
        return 1

    # Build the plotter
    plotter = onc_parser.make_plotter()

    # Table of methods
    callables = collections.OrderedDict([
        ("wp", plotter.plot_waves_and_projs),
        ("dp", plotter.plot_dens_and_pots),
        ("lc", plotter.plot_atanlogder_econv),
        ("df", plotter.plot_den_formfact),
    ])

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
    out_path = options.filename
    onc_parser = OncvOutputParser(out_path)
    onc_parser.scan()
    if not onc_parser.run_completed:
        cprint("oncvpsp output is not complete. Exiting", "red")
        return 1

    # Generate json files with oncvpsp results.
    print(json.dumps(onc_parser.to_dict, indent=-1))
    return 0


def oncv_run(options):
    """Run oncvpsp, generate djrepo file, plot results. Requires input file."""

    calc_type = dict(nor= "non-relativistic",
                     sr = "scalar-relativistic",
                     fr="fully-relativistic")[options.rel]

    from pseudo_dojo.ppcodes import OncvGenerator
    in_path = options.filename
    oncv_ppgen = OncvGenerator.from_file(in_path, calc_type, workdir=None)
    #print(oncv_ppgen)
    #oncv_ppgen.start()
    #oncv_ppgen.wait()

    #if retcode != 0
    #   cprint("oncvpsp returned %s. Exiting" % retcode, "red")
    #   return retcode

    out_path = in_path.replace(".in", ".out")
    """
    onc_parser = OncvOutputParser(out_path)
    onc_parser.scan()
    if not onc_parser.run_completed:
        cprint("oncvpsp output is not complete. Exiting", "red")
        return 1

    if options.psp8:
        # Generate psp8 and djson files.
        root, _ = os.path.splitext(out_path)
        psp8_path = root + ".psp8"
        djson_path = root + ".djson"

        # Extract psp8 files from the oncvpsp output and write it to file.
        s = onc_parser.get_pseudo_str(devel=True)
        with open(psp8_path, "wt") as fh:
            fh.write(s)

        # Try to read pseudo from the files just generated.
        pseudo = dojopseudo_from_file(psp8_path)
        if pseudo is None:
            cprint("Cannot parse psp8 files: %s" % psp8_path, "red")
            return 1

        # Write djson file.
        #report = DojoReport.new_pseudo(pseudo)
        #report.json_write(djson_path)
        return 0
    """


def main():

    def str_examples():
        return """\
Usage example:
    oncv.py run H.in                  ==> Run oncvpsp.
    oncv.py plot H.out                ==> Plot oncvpsp generation results for pseudo H.psp8
    oncv.py json H.out                ==> Generate JSON file.
"""

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg: sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    # Parent parser implementing commong options
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
    p_run.add_argument("--rel", default="sr", help="Relativistic treatment: nor, sr, fr")

    # Create the parsers for the sub-commands
    p_plot = subparsers.add_parser('plot', parents=[copts_parser], help=oncv_plot.__doc__)
    p_plot.add_argument("-p", "--plot-mode", default="slide",
                        help=("Quantity to plot. Possible values: %s" %
                              str(["slide", "wp, dp, lc"] + PseudoGenDataPlotter.all_keys) + "\n"
                              "wp --> wavefunctions and projectors\n" +
                              "dp --> densities and potentials\n" +
                              "lc --> atan(logder) and convergence wrt ecut\n" +
                              "df --> density form factor"))

    p_json = subparsers.add_parser('json', parents=[copts_parser], help=oncv_json.__doc__)

    #parser.add_argument("-8", "--psp8", action="store_true", default=False,
    #                    help="produce a .psp8 file with initial dojo report and exit")
    #parser.add_argument("-d", "--devel", action="store_true", default=False,
    #                    help="put only two energies in the ecuts list for testing for developing the pseudo")

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
