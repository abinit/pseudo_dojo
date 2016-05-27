#!/usr/bin/env python
"""Script to analyze and plot the results produced by ONCVPSP code."""
from __future__ import division, print_function, unicode_literals

import sys
import collections
import argparse

from monty.termcolor import cprint 
from pseudo_dojo.core.pseudos import dojopseudo_from_file
from pseudo_dojo.ppcodes.oncvpsp import OncvOutputParser, PseudoGenDataPlotter


def main():
    parser = argparse.ArgumentParser() #formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('filename', default="", help="Path to the output file")

    parser.add_argument("-p", "--plot-mode", default="slide",
                        help=("Quantity to plot. Possible values: %s" %
                              str(["slide", "wp, dp, lc"] + PseudoGenDataPlotter.all_keys) + "\n"
                              "wp --> wavefunctions and projectors\n" +
                              "dp --> densities and potentials\n" +
                              "lc --> atan(logder) and convergence wrt ecut\n" + 
                              "df --> density form factor"))

    parser.add_argument("-j", "--json", action="store_true", default=False, 
                        help="Produce a string with the results in a JSON dictionary and exit")

    parser.add_argument("-8", "--psp8", action="store_true", default=False, 
                        help="produce a .psp8 file with initial dojo report and exit")

    parser.add_argument("-d", "--devel", action="store_true", default=False,
                        help="put only two energies in the ecuts list for testing for develloping the pseudo")

    options = parser.parse_args()

    onc_parser = OncvOutputParser(options.filename)
    onc_parser.scan()
    if not oncv_parser.completed:
        cprint("oncvpsp output is not complete. Exiting", "red")
        return 1

    if options.psp8:
        # Generate psp8 and djson files.
        root, _ = os.path.splitext(options.filename)
        psp8_path = root + ".psp8"
        djson_path = root + ".djson"

        # Extract psp8 files from out and write it to file.
        s = onc_parser.get_pseudo_str(options.devel=True)
        with open(psp8_path, "wt") as fh:
            fh.write(s)

        #report = DojoReport.new_from_file()
        #report.json_write(djson_path)

        # Try to read pseudo from the files just generated.
        pseudo = dojopseudo_from_file(psp8_path)
        if pseudo is None:
            cprint(("Cannot parse psp8 files: %s" % psp8_path, "red")
            return 1

        return 0      

    if options.json:
        # Generate json files with oncvpsp results.
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
        ("df", plotter.plot_den_formfact),
    ])

    #plotter.plot_radial_wfs()
    #plotter.plot_projectors()
    #plotter.plot_potentials()
    #plotter.plot_der_potentials()
    #for order in [1,2,3,4]:
    #    plotter.plot_der_densities(order=order)
    #plotter.plot_densities()
    #plotter.plot_densities(timesr2=True)
    #plotter.plot_den_formfact()
    #return 0

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
