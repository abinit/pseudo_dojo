#!/usr/bin/env python
"""Script to check pseudopotential files."""
from __future__ import division, print_function, unicode_literals

import sys
import argparse

__author__ = "Matteo Giantomassi"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"

def check(options):
    """Check pseudopotential"""
    return 0

#def regenerate(options):
#    return 0

def main():
    def str_examples():
        return """\n
        Example usage:\n
            pscheck.py check       => check the pseudopotential.
            pscheck.py regenerate  => Generate new md5 hash values.
        """

    def show_examples_and_exit(err_msg=None, error_code=0):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg:
            sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    parser = argparse.ArgumentParser(epilog=str_examples())

    parser.add_argument('-V', '--version', action='version', version="%(prog)s version " + __version__)
    parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
                         help='verbose, can be supplied multiple times to increase verbosity')

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help')

    # Subparser for check command.
    p_check = subparsers.add_parser('check', help='check pseudos.')

    # Subparser for regenerate command.
    #p_reg = subparsers.add_parser('regenerate', help='Regenerate new checksums.')

    # Parse command line.
    try:
        options = parser.parse_args()
    except Exception as exc: 
        show_examples_and_exit(error_code=1)

    # Dispatch
    return globals()[options.command](options)

    return 0

if __name__ == "__main__":
    try:
        do_prof = sys.argv[1] == "prof"
        if do_prof: sys.argv.pop(1)
    except Exception:
        do_prof = False

    if do_prof:
        import pstats, cProfile
        cProfile.runctx("main()", globals(), locals(), "Profile.prof")
        s = pstats.Stats("Profile.prof")
        s.strip_dirs().sort_stats("time").print_stats()
    else:
        sys.exit(main())
