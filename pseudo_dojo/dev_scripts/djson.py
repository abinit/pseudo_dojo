#!/usr/bin/env python
"""
Script to validate/regenerate djson files.
"""
from __future__ import division, print_function, unicode_literals

import sys
import argparse
import json

from pprint import pprint
from pseudo_dojo.core.pseudos import DojoTable

__author__ = "Matteo Giantomassi"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"

def djson_new(options):
    """Create new djson file. Print document to stdout."""
    # Build full table
    table = DojoTable.from_dojodir(options.top)
    djson = table.to_djson()
    print(json.dumps(djson, indent=-1))
    return 0


def djson_validate(options):
    """Validate djson file."""
    table = DojoTable.from_djson_file(options.djson_path)
    print(table)
    errors = table.validate()
    if errors:
        pprint(errors)

    return len(errors)


def main():
    def str_examples():
        return """\
Usage example: 
    djson.py  new [DIR]        => Generate new djson file.
    djson.py  validate djson   => Validate djson file.
"""

    def show_examples_and_exit(err_msg=None, error_code=0):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg:
            sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=str_examples(), formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-V', '--version', action='version', version="%(prog)s version " + __version__)
    parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
                         help='verbose, can be supplied multiple times to increase verbosity')

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help')

    # Subparser for new command.
    p_new = subparsers.add_parser('new', help='Generate new djson file.')
    p_new.add_argument("top", nargs="?", default=".", help="Directory with pseudos")

    # Subparser for validate command.
    p_validate = subparsers.add_parser('validate', help='Validate the djson file.')
    p_validate.add_argument("djson_path", help="dsjon file")

    # Parse command line.
    try:
        options = parser.parse_args()
    except Exception as exc: 
        show_examples_and_exit(error_code=1)

    # Dispatch
    return globals()["djson_" + options.command](options)


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
