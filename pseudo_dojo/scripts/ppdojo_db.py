#!/usr/bin/env python
"""Simple search engine for the pseudodojo databases."""
from __future__ import division, print_function

import sys
import numpy as np

from argparse import ArgumentParser

from pseudo_dojo import pseudodojo_database

__author__ = "Matteo Giantomassi"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"

##########################################################################################

def str_examples():
    examples = """Example usage:\n
      ppdojo_db.py nc_find Si  => Show all norm-conserving pseudos for Si.
    """
    return examples


def show_examples_and_exit(err_msg=None, error_code=0):
    "Display the usage of the script."
    sys.stderr.write(str_examples())
    if err_msg: 
        sys.stderr.write("Fatal Error\n" + err_msg + "\n")
    sys.exit(error_code)

##########################################################################################

def main():
    parser = ArgumentParser(epilog=str_examples())

    #parser.add_argument('--version', action='version', version="ppdojo_db.py " + __version__)

    parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
                         help='verbose, can be supplied multiple times to increase verbosity')

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help')

    # Subparser for show command.
    p_show = subparsers.add_parser('show', help='Show info on the database.')

    # Subparser for nc_find command.
    p_find = subparsers.add_parser('nc_find', help='Find pseudopotentials in the database.')

    p_find.add_argument('-xc', '--xc-type', default="GGA", help="XC functional type. (default GGA).")  

    p_find.add_argument('-t', '--table-name', default=None, help="Table name (default None).")  

    p_find.add_argument('-s', '--sort', default=None, help="Sort pseudos according to this attribute.")  

    p_find.add_argument('symbol', nargs='?', help='Chemical symbol')


    ###################################
    # Parse command line and dispatch #
    ###################################
    options = parser.parse_args()

    verbose = options.verbose

    #reload_pseudodojo_database():
    ppdb = pseudodojo_database()

    if options.command == "show":
        #pp_type = options.pp_type
        #xc_type = options.xc_type
        pp_type = "NC"
        xc_type = "GGA"
        ppdb.show(pp_type, xc_type, verbose=verbose)

    if options.command == "nc_find":
        symbol = options.symbol
        xc_type = options.xc_type

        pseudos = ppdb.nc_pseudos(symbol, xc_type, table_name=options.table_name)

        # Sort pseudos.
        if options.sort:
            attrs = []
            for i, p in pseudos:
                try:
                    a = getattr(p, options.sort)
                except AttributeError:
                    a = np.inf
                attrs.append((i,a))

            # Sort attrs, then shuffle pseudos.
            attrs = sorted(attrs, key=lambda t:t[1])
            pseudos = [pseudos[a[0]] for a in attrs]

        for p in pseudos:
            print()
            print(p)
            if p.has_dojo_report: 
                print(p.dojo_report)
            print()

    return 0

################################################################################

if __name__ == "__main__":
    sys.exit(main())
