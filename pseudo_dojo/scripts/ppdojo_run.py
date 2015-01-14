#!/usr/bin/env python
"""This scripts drives the execution of the pseudo dojo tests."""
from __future__ import division, print_function, unicode_literals

import sys
import argparse

from pprint import pprint
from pymatgen.io.abinitio import TaskManager
from pymatgen.io.abinitio.pseudos import PseudoTable 
from pseudo_dojo.dojo import Dojo, HintsAndGbrvDojo

__author__ = "Matteo Giantomassi"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"


def main():

    def str_examples():
        examples = """
Usage Example:\n
    ppdojo_run.py build Si.fhi  => Build pseudo_dojo flow for Si.fhi
\n"""
        return examples

    def show_examples_and_exit(error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        sys.exit(error_code)

    parser = argparse.ArgumentParser(epilog=str_examples())

    #parser.add_argument('-l', '--max-level', type=int, default=0,  help="Maximum DOJO level (default 0 i.e. ecut hints).")

    parser.add_argument('-m', '--manager', type=str, default=None,  help="Manager file")

    parser.add_argument('--loglevel', default="ERROR", type=str,
                        help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")

    parser.add_argument('pseudos', nargs='+', help='List of pseudopotential files.')

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    # Subparser for single command.
    p_build = subparsers.add_parser('build', help="Build dojo.")

    # Subparser for single command.
    p_report = subparsers.add_parser('report', help="Show DOJO_REPORT.")

    try:
        options = parser.parse_args()
    except:
        show_examples_and_exit(1)

    # loglevel is bound to the string value obtained from the command line argument. 
    # Convert to upper case to allow the user to specify --loglevel=DEBUG or --loglevel=debug
    import logging
    numeric_level = getattr(logging, options.loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % options.loglevel)
    logging.basicConfig(level=numeric_level)

    pseudos = PseudoTable(options.pseudos) 
    if options.manager is None:
        manager = TaskManager.from_user_config()
    else:
        manager = TaskManager.from_file(options.manager)

    if options.command == "build":
        dojo = HintsAndGbrvDojo(manager=manager)
        #dojo = Dojo(manager=manager)

        for pseudo in pseudos:
            dojo.add_pseudo(pseudo)

        #dojo.build()
        dojo.start()

    elif options.command == "report":
        for pseudo in pseudos:
            report = pseudo.report
            report.print_table()

    return 0

if __name__ == "__main__":
    sys.exit(main())
