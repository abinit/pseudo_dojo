#!/usr/bin/env python
"""Script used by maintainers to check djrepo files."""
from __future__ import division, print_function, unicode_literals

import sys
import os
import argparse
import glob

from monty.functools import prof_main
from monty.termcolor import cprint 
from monty.os.path import find_exts
from pseudo_dojo.core.pseudos import dojopseudo_from_file


def _fix_djrepo(pp_filepath):
    """
    This is a maintentance tool:
        #. Regenerate the md5 value in the DojoReport file.
        #. Replace norm-conserving with NC. 
    """
    pseudo = dojopseudo_from_file(pp_filepath)
    if pseudo is None:
        print("Error while parsing %s" % pp_filepath)
        return 1

    # Change md5
    pseudo.dojo_report["md5"] = pseudo.compute_md5()
    if pseudo.dojo_report["pseudo_type"] == "norm-conserving":
        pseudo.dojo_report["pseudo_type"] = "NC"

    pseudo.dojo_report.json_write(pseudo.djrepo_path)
    return 0


def djrepo_check(options):
    """Check djrepo files."""
    return 0

#def djrepo_regmd5(options):
#    return 0


#def djrepo_convert(options):
#    new_version = "2.0"


@prof_main
def main():
    def str_examples():
        return """\
Example usage:
    djrepo.py check     => Check djrepo files.
    djrepo.py convert   => Convert djrepo to new format.
    djrepo.py recalc    => Recalculate deltafactor in djrepo
    djrepo.py md5reg    => Generate new md5 hash values.
"""

    def show_examples_and_exit(err_msg=None, error_code=0):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg:
            sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    # Parent parser for commands that need to know on which subset of pseudos we have to operate.
    copts_parser = argparse.ArgumentParser(add_help=False)
    copts_parser.add_argument('paths', nargs="+", help="Pseudopotential file or directory containing pseudos")
    copts_parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
                         help='verbose, can be supplied multiple times to increase verbosity')

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=str_examples(), formatter_class=argparse.RawDescriptionHelpFormatter)

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help')

    # Subparser for check command.
    p_check = subparsers.add_parser('check', parents=[copts_parser], help=djrepo_check.__doc__)

    # Subparser for regenerate command.
    #p_reg = subparsers.add_parser('md5reg', parents=[copts_parser], help=md5reg.__doc__)

    # Subparser for convert command.
    #p_reg = subparsers.add_parser('convert', parents=[copts_parser], help=md5reg.__doc__)

    # Parse command line.
    try:
        options = parser.parse_args()
    except Exception as exc: 
        show_examples_and_exit(error_code=1)

    # Get paths
    def get_paths(options):
        """
        Find and retur list of djson paths from options
        Accepts filepaths or directory.
        """
        paths = options.paths

        if len(paths) == 1:
            # Handle directory argument
            if os.path.isdir(paths[0]):
                top = os.path.abspath(paths[0])
                paths = find_exts(top, ["djrepo"], exclude_dirs="_*")
            # Handle glob syntax e.g. "./*.psp8" 
            elif "*" in paths[0]:
                paths = [f for f in glob.glob(paths[0]) if f.endswith(".djrepo")]
        else:
            paths = [p for p in paths if p.endswith(".djrepo")]

        return paths

    options.paths = get_paths(options)
    if not options.paths:
        cprint("Empty djrepo list. Returning", "magenta")
        return 1

    if options.verbose:
        print("Got %s djrepo files" % len(options.paths))
        if options.verbose > 1:
            for i, p in enumerate(options.paths): print("[%d] %s" % (i, p))

    # This is to regenerate the md5 files in the djrepo files.
    """
    from pseudo_dojo.pseudos import _fix_djrepo
    for pseudo in options.pseudos:
        _fix_djrepo(pseudo.filepath)
    return 0
    """

    # Dispatch
    return globals()["djrepo_" + options.command](options)


if __name__ == "__main__":
    sys.exit(main())
