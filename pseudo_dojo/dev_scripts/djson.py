#!/usr/bin/env python
"""Script to validate/regenerate djson files."""
from __future__ import division, print_function, unicode_literals

import sys
import argparse
import json

from monty.functools import prof_main
from monty.termcolor import cprint 
from pseudo_dojo.core.pseudos import DojoTable, OfficialDojoTable


def djson_new(options):
    """Create new djson file. Print document to stdout."""
    # Build full table
    table = DojoTable.from_dojodir(options.top)
    djson = table.to_djson()
    print(json.dumps(djson, indent=-1))
    return 0


def djson_validate(options):
    """Validate djson file."""
    table = OfficialDojoTable.from_djson_file(options.djson_path)
    #print(table)

    md5dict = {p.basename: p.md5 for p in table}
    errors = table.dojo_find_errors(md5dict, require_hints=False)

    if errors:
        cprint("dojo_find_errors returned %s errors" % len(errors), "red")
        if not options.verbose:
            print("Use --verbose to show errors")
        else:
            for i, e in enumerate(errors): 
                print("[%s]" % i, e)

    return len(errors)


@prof_main
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

    #parser.add_argument('-V', '--version', action='version', version="%(prog)s version " + __version__)
    parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
                         help='verbose, can be supplied multiple times to increase verbosity')

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help')

    # Subparser for new command.
    p_new = subparsers.add_parser('new', help=djson_new.__doc__)
    p_new.add_argument("top", nargs="?", default=".", help="Directory with pseudos")

    # Subparser for validate command.
    p_validate = subparsers.add_parser('validate', help=djson_validate.__doc__)
    p_validate.add_argument("djson_path", help="djson file")

    # Parse command line.
    try:
        options = parser.parse_args()
    except Exception as exc: 
        show_examples_and_exit(error_code=1)

    # Dispatch
    return globals()["djson_" + options.command](options)

if __name__ == "__main__":
    sys.exit(main())
