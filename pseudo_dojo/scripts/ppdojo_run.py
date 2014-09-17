#!/usr/bin/env python
"""This scripts drives the execution of the pseudo dojo tests."""
from __future__ import division, print_function, unicode_literals

import sys
import argparse

from pprint import pprint
from pymatgen.io.abinitio import TaskManager
from pseudo_dojo.dojo import Dojo, HintsAndGbrvDojo, DojoReport

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

    # Decorate argparse classes to add portable support for aliases in add_subparsers
    class MyArgumentParser(argparse.ArgumentParser):
        def add_subparsers(self, **kwargs):
            new = super(MyArgumentParser, self).add_subparsers(**kwargs)
            # Use my class
            new.__class__ = MySubParserAction
            return new
                                                                                                                    
    class MySubParserAction(argparse._SubParsersAction):
        def add_parser(self, name, **kwargs):
            """Allows one to pass the aliases option even if this version of ArgumentParser does not support it."""
            try:
                return super(MySubParserAction, self).add_parser(name, **kwargs)
            except Exception as exc:
                if "aliases" in kwargs: 
                    # Remove aliases and try again.
                    kwargs.pop("aliases")
                    return super(MySubParserAction, self).add_parser(name, **kwargs)
                else:
                    # Wrong call.
                    raise exc

    parser = MyArgumentParser(epilog=str_examples())

    #parser.add_argument('-l', '--max-level', type=int, default=0, 
    #                    help="Maximum DOJO level (default 0 i.e. ecut hints).")

    parser.add_argument('--loglevel', default="ERROR", type=str,
                        help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")

    parser.add_argument('pseudos', nargs='+', help='List of pseudopotential files.')

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    # Subparser for single command.
    p_build = subparsers.add_parser('build', aliases=["b"], help="Build dojo.")

    # Subparser for single command.
    p_report = subparsers.add_parser('report', aliases=["r"], help="Show DOJO_REPORT.")

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

    pseudos = options.pseudos
    manager = TaskManager.from_user_config()

    if options.command == "build":
        #dojo = Dojo(manager=manager)
        dojo = HintsAndGbrvDojo(manager=manager)

        for pseudo in pseudos:
            dojo.add_pseudo(pseudo)

        dojo.build()

    elif options.command == "report":
        for pseudo in pseudos:
            report = DojoReport.from_file(pseudo)
            report.print_table()

    return 0

if __name__ == "__main__":
    sys.exit(main())
