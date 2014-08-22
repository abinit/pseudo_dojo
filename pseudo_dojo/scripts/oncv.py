#!/usr/bin/env python
"""Script to analyze and plot the results produced by ONCVPSP code."""
from __future__ import print_function, division

import sys
import collections

from pseudo_dojo.ppcodes.oncvpsp import  OncvOuptputParser, PseudoGenDataPlotter


def main():
    import argparse
    parser = argparse.ArgumentParser() #formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('filename', default="", help="Path to the output file")

    parser.add_argument("-p", "--plot-mode", default="slide",
                        help=("Quantity to plot. Possible values: %s" %
                              str(["slide", "wp, dp, lc"] + PseudoGenDataPlotter.all_keys) + "\n"
                              "wp --> wavefunctions and projectors\n" +
                              "dp --> densities and potentials\n" +
                              "lc --> atan(logder) and convergence wrt ecut"))

    options = parser.parse_args()

    onc_parser = OncvOuptputParser(options.filename)
    onc_parser.scan()
    print(onc_parser)

    # Build the plotter
    plotter = onc_parser.make_plotter()

    # Table of methods
    callables = collections.OrderedDict([
        ("wp", plotter.plot_waves_and_projs),
        ("dp", plotter.plot_dens_and_pots),
        ("lc", plotter.plot_atanlogder_econv),
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

    return 0


if __name__ == "__main__":
    sys.exit(main())

