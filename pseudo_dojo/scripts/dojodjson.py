#!/usr/bin/env python
"""Script to validate/regenerate djson files."""
import sys
import os
import argparse
import json

from monty.functools import prof_main
from monty.termcolor import cprint
from pseudo_dojo.core.pseudos import DojoTable, OfficialDojoTable


def djson_new(options, stream=sys.stdout):
    """Create a new djson file from directory or from text file. Print document to stdout."""
    # Build full table (either from directory or text file.
    if not os.path.isdir(options.top):
        # Assume top is a text file
        table = DojoTable.from_txtfile(options.top)
    else:
        table = DojoTable.from_dojodir(options.top, exclude_wildcard=options.exclude)
    djson = table.to_djson(options.verbose, ignore_dup=False)
    print(json.dumps(djson, indent=1), file=stream)
    return 0


def djson_validate(options):
    """Validate djson file."""
    table = OfficialDojoTable.from_djson_file(options.djson_path)
    if options.verbose > 1: print(table)

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
    djson.py  new [DIR]        => Generate new djson file from pseudos inside directory.
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

    # Parent parser for common options.
    copts_parser = argparse.ArgumentParser(add_help=False)
    copts_parser.add_argument('--loglevel', default="ERROR", type=str,
                               help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")
    copts_parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
                               help='verbose, can be supplied multiple times to increase verbosity')

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help')

    # Subparser for new command.
    p_new = subparsers.add_parser('new', parents=[copts_parser], help=djson_new.__doc__)
    p_new.add_argument("--exclude", default=None, type=str,
                       help=("Exclude files mathing these pattern. Example"
                             "exclude=\"*_r.psp8|*.xml\" selects only those files that do not end with _r.psp8 or .xml"
                      ))
    p_new.add_argument("top", nargs="?", default=".", help="Directory with pseudos or text file")

    # Subparser for validate command.
    p_validate = subparsers.add_parser('validate', parents=[copts_parser], help=djson_validate.__doc__)
    p_validate.add_argument("djson_path", help="djson file")

    # Parse command line.
    try:
        options = parser.parse_args()
    except Exception as exc:
        show_examples_and_exit(error_code=1)

    # loglevel is bound to the string value obtained from the command line argument.
    # Convert to upper case to allow the user to specify --loglevel=DEBUG or --loglevel=debug
    import logging
    numeric_level = getattr(logging, options.loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % options.loglevel)
    logging.basicConfig(level=numeric_level)

    # Dispatch
    return globals()["djson_" + options.command](options)


if __name__ == "__main__":
    sys.exit(main())
