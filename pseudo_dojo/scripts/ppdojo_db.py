#!/usr/bin/env python
"""Search engine for the pseudodojo databases."""
from __future__ import division, print_function

import sys

from argparse import ArgumentParser

from pseudo_dojo import pseudodojo_database

__author__ = "Matteo Giantomassi"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"

##########################################################################################
# Helper functions.


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
    parser = ArgumentParser()

    #parser.add_argument('--version', action='version', version="ppdojo_db.py " + __version__)

    parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
                         help='verbose, can be supplied multiple times to increase verbosity')

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help')

    # Subparser for show command.
    p_show = subparsers.add_parser('show', help='Show info on the database.')

    # Subparser for show command.
    p_find = subparsers.add_parser('nc_find', help='Find pseudopotentials in the database.')

    p_find.add_argument('-xc', '--xc-type', default="GGA", help="XC functional type. (default GGA).")  

    p_find.add_argument('-t', '--table-name', default=None, help="Table name (default None).")  

    p_find.add_argument('symbol', nargs='?', help='Chemical symbol')


    ###################################
    # Parse command line and dispatch #
    ###################################
    options = parser.parse_args()

    verbose = options.verbose

    #reload_pseudodojo_database():
    ppdb = pseudodojo_database()

    if options.command == "show":
        ppdb.show(verbose)

    if options.command == "nc_find":
        symbol = options.symbol
        xc_type = options.xc_type

        pseudos = ppdb.nc_pseudos(symbol, xc_type, table_name=options.table_name)

        for p in pseudos:
            print()
            print(p)
            print()

    return 0

################################################################################

if __name__ == "__main__":
    sys.exit(main())
