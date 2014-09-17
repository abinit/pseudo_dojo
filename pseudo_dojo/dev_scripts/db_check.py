#!/usr/bin/env python
"""Script to validate/regenerate the database of md5 values."""
from __future__ import division, print_function, unicode_literals

import sys

from argparse import ArgumentParser
from pseudo_dojo.pseudos.database import validate_checksums, regenerate_checksums

__author__ = "Matteo Giantomassi"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"

##########################################################################################
# Helper functions.


def str_examples():
    examples = """Example usage:\n
      db_check.py  validate    => Validate the pseudopotentials.
      db_check.py  regenerate  => Generate new md5 hash values.
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

    # Subparser for validate command.
    p_validate = subparsers.add_parser('validate', help='Validate the pseudopotential database.')

    # Subparser for regenerate command.
    p_reg = subparsers.add_parser('regenerate', help='Regenerate new checksums.')

    ###################################
    # Parse command line and dispatch #
    ###################################
    options = parser.parse_args()

    verbose = options.verbose

    if options.command == "validate":
        isok = validate_checksums(verbose=verbose)
        print("isok %s " % isok)

    if options.command == "regenerate":
        ans = raw_input("This command will regenerate the md5 hash values. Enter Y to continue... ")
        if ans == "Y":
            regenerate_checksums(verbose=verbose)

    return 0

################################################################################

if __name__ == "__main__":
    sys.exit(main())
