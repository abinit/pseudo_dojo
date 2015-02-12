#!/usr/bin/env python
"""Script to analyze and plot the results produced by ONCVPSP code."""
from __future__ import division, print_function, unicode_literals

import sys
import collections
import argparse

from pseudo_dojo.ppcodes.oncvpsp import OncvOutputParser, PseudoGenDataPlotter


def main():
    parser = argparse.ArgumentParser() #formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('filename', default="", help="Path to the output file")

    parser.add_argument("-p", "--plot-mode", default="slide",
                        help=("Quantity to plot. Possible values: %s" %
                              str(["slide", "wp, dp, lc"] + PseudoGenDataPlotter.all_keys) + "\n"
                              "wp --> wavefunctions and projectors\n" +
                              "dp --> densities and potentials\n" +
                              "lc --> atan(logder) and convergence wrt ecut"))

    parser.add_argument("-j", "--json", action="store_true", default=False, 
                        help="Produce a string with the results in a JSON dictionary and exit")

    options = parser.parse_args()

    onc_parser = OncvOutputParser(options.filename)
    onc_parser.scan()

    if options.json:
        import json
        print(json.dumps(onc_parser.to_dict, indent=4))
        return 0

    # Build the plotter
    plotter = onc_parser.make_plotter()

    # Table of methods
    callables = collections.OrderedDict([
        ("wp", plotter.plot_waves_and_projs),
        ("dp", plotter.plot_dens_and_pots),
        ("lc", plotter.plot_atanlogder_econv),
    ])

    #plotter.plot_radial_wfs()
    #plotter.plot_projectors()
    #plotter.plot_potentials()
    #plotter.plot_der_potentials()
    #plotter.plot_densities()
    #plotter.plot_der_densities(order=1)
    #plotter.plot_der_densities(order=2)
    plotter.plot_der_densities(order=4)
    return

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

