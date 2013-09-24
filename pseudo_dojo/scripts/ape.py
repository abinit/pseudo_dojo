#!/usr/bin/env python
from __future__ import division, print_function

import os
import sys
import glob
import argparse

from pseudo_dojo.refdata.nist import nist_database 
from pseudo_dojo.ppcodes.ape import *

__version__ = "0.1"

##########################################################################################

def str_examples():
    examples = """Example usage:\n
      ape.py nist Si Ca   => Show NIST LDA data for Si and Ca (http://physics.nist.gov/PhysRefData/DFTdata/)
    """
    return examples


def show_examples_and_exit(err_msg=None, error_code=0):
    """Display the usage of the script."""
    sys.stderr.write(str_examples())
    if err_msg: 
        sys.stderr.write("Fatal Error\n" + err_msg + "\n")
    sys.exit(error_code)


def issymbol(string):
    """True if string is a known element symbol."""
    return string in nist_database.allsymbols


def show_nistdata(options):
    """Handle nist command."""
    for symbol in options.symbols:
        if symbol.endswith("+"):
            iontype = "Cation"
            raise NotImplementedError("Cations not supported")
        else:
            iontype = "Neutral"
            entry = nist_database.get_neutral_entry(symbol)
                                                               
        print("")
        print(entry)
        print("")


def plot(options):
    "Handle plot command"
    dirpath = options.dirpath

    if "w" in options.plots:
        ape_plot_waves(dirpath, savefig=None)
                                                
    if "l" in options.plots:
        ape_plot_logders(dirpath, savefig=None)

    if "p" in options.plots:
        ape_plot_potentials(dirpath, savefig=None)

    if "d" in options.plots:
        ape_plot_densities(dirpath, savefig=None)


def analyze_llocal(options):
    """Handle llocal command."""
    # Create an InputGenerator to facilitate the modification of the input file.
    inpgen = ApeInputGenerator.from_template(options.template)

    # Init the AE solver from the template
    ae_solver = ApeAeSolver(workdir, inpgen, verbose=verbose)
                                                                        
    if verbose: 
        ae_solver.show_input()
                               
    if not dry_run:
        ae_solver.solve(remove_wd=remove_wd)

        # Define the parameters for the pseudization.
        # For each possible local L:
        #    1) Pseudize.
        #    2) Check ghosts
        #    3) Plot wfs and logders

        pp_generators = []
        #for llocal in range(-1, 5, 1):
        for llocal in range(0, 1, 1):
            inpgen.reset()
            inpgen.set_llocal(llocal)
            #inpgen.set_ppcomponents()
            #inpgen.set_correction()

            #pp_components = ApePPComponents.from_strings("3s|1.2|tm", "3p|1.27|tm", "3d|1.5|tm", "4f|1.9|tm")
            #pp_setup = ApePPSetup(pp_components, core_correction=0, llocal=llocal)

            pp_workdir = os.path.join(workdir, "ppgen_loc%d" % llocal)

            pp_gen = ApePseudoGenerator(pp_workdir, inpgen, ae_solver, verbose=verbose)

            pp_generators.append(pp_gen)

        for pp_gen in pp_generators:
            pp_gen.pseudize(remove_wd=remove_wd)
            #if pp_gen.ghosts:
            #    print("Detected ghosts for states %s" % pp_gen.ghosts.keys())


##########################################################################################


def main():
    parser = argparse.ArgumentParser(epilog=str_examples(),
                                     formatter_class=argparse.RawDescriptionHelpFormatter,)

    parser.add_argument('--version', action='version', version="ape.py " + __version__)

    parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
                        help='verbose, can be supplied multiple times to increase verbosity')  

    parser.add_argument('-d', '--dry-run', action='store_true', default=False, help='Dry run mode')  

    parser.add_argument('-r', '--remove_wd', action='store_true', default=False, 
                        help="Remove working dirs if they already exist")  

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help')

    # AE subparser
    #p_ae = subparsers.add_parser('ae', help='Solve the all-electron problem')

    #p_ae.add_argument('aconf_str', metavar='STRING', default= "", help="String with the all-electron configuration" + 
    #                  "It can be either element symbol (for neutral configuration) or string in the form [He] 2s2...") 

    # Automatic generation of a (very) simple input file.
    #p_autogen = subparsers.add_parser('autogen', help='Solve the all-electron problem')

    #p_autogen.add_argument('aconf_str', metavar='STRING', default= "", help="String with the all-electron configuration" + 
    #                      "It can be either element symbol (for neutral configuration) or string in the form [He] 2s2...") 

    # Subparser for Pseudization with template.
    p_pst = subparsers.add_parser('pst', help='Pseudization with template input file')

    p_pst.add_argument('template', metavar='STRING', default= "", help="Path to the template file")

    # Plot subparser
    p_plot = subparsers.add_parser('plot', help='Plot data with matplotlib')

    p_plot.add_argument('-p', '--plots', metavar='STRING', default= "w", help="Quantities to plot: " + 
                        "w for wavefunctions, l for logarithmic derivatives, p for potentials, d for densities (DEFAULT: w)")

    p_plot.add_argument('dirpath', metavar='DIRPATH', default = ".", help="Path to the directory containing APE output file") 

    # Analysis of LLocal option.
    #p_llocal = subparsers.add_parser('llocal', help='Analyze the choice of the angular momentum for the local part')

    #p_llocal.add_argument('template', metavar='STRING', default= "", help="Path to the template file")

    # Access the NIST database
    p_nist = subparsers.add_parser('nist', help='Retrieve AE results from the NIST database')

    p_nist.add_argument('symbols', nargs="+", help="List of element symbols")

    ###################################
    # Parse command line and dispatch #
    ###################################
    options = parser.parse_args()
                                  
    verbose = options.verbose
    dry_run  = options.dry_run
    remove_wd = options.remove_wd
                                  


    #if options.command == "autogen":
    #    aconf_str = options.aconf_str
    #                                                                      
    #    if issymbol(aconf_str):
    #        # AE configuration: neutral silicon + empty 3d-4f states
    #        aconf = ApeAtomicConfiguration.neutral_from_symbol(aconf_str)
    #    else:
    #        raise NotImplementedError("")

    #    #inpgen = ApeInputGenerator.from_objects()
    #    #print(inpget.get_strlist())
                                                                          
    #if options.command == "ae":
    #    raise NotImplementedError("this part must be re-tested")
    #    aconf_str = options.aconf_str

    #    if issymbol(aconf_str):
    #        # AE configuration: neutral silicon + empty 3d-4f states
    #        aconf = ApeAtomicConfiguration.neutral_from_symbol(aconf_str)

    #    else:
    #        raise NotImplementedError("")
    #        aconf.add_state(n=3, l="d", occ=0.0)
    #        aconf.add_state(n=4, l="f", occ=0.0)

    #    # Solve AE problem.
    #    ae_solver = ApeAeSolver(workdir, aconf, verbose=verbose)

    #    if verbose: ae_solver.show_input()

    #    if not dry_run:
    #        ae_solver.solve(remove_wd=remove_wd)

    if options.command == "pst":

        # Create an InputGenerator to facilitate the modification of the input file.
        inpgen = ApeInputGenerator.from_template(options.template)

        # Construct the name of the working directory. Ex: APERUN_1, APERUN_2 ...
        prefix = "APERUN_"
        nums = [int(dirname[len(prefix):]) for dirname in glob.glob(prefix + "*")]
        max_num = max(nums) if nums else 0
        workdir = prefix + str(max_num+1)

        # Init the AE solver from the template
        ae_solver = ApeAeSolver(workdir, inpgen, verbose=verbose)

        if not dry_run:
            ae_solver.solve(remove_wd=remove_wd)

            # Define the parameters for the pseudization.
            # For each possible local L:
            #    1) Pseudize.
            #    2) Check ghosts
            #    3) Plot wfs and logders

            pp_generators = []
            #for llocal in range(-1, 5, 1):
            for llocal in range(0, 1, 1):

                inpgen.reset()
                inpgen.set_llocal(llocal)
                #inpgen.set_ppcomponents()
                #inpgen.set_correction()

                #pp_components = ApePPComponents.from_strings("3s|1.2|tm", "3p|1.27|tm", "3d|1.5|tm", "4f|1.9|tm")
                #pp_setup = ApePPSetup(pp_components, core_correction=0, llocal=llocal)

                pp_workdir = os.path.join(workdir, "ppgen_loc%d" % llocal)

                pp_gen = ApePseudoGenerator(pp_workdir, inpgen, ae_solver, verbose=verbose)

                pp_generators.append(pp_gen)

            for pp_gen in pp_generators:
                pp_gen.pseudize(remove_wd=remove_wd)


    if options.command == "llocal":
        analyze_llocal(options)

    if options.command == "plot":
        plot(options)

    if options.command == "nist":
        show_nistdata(options)

    return 0

##########################################################################################

if __name__ == "__main__":
    sys.exit(main())
