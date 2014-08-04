#!/usr/bin/env python
"""This scripts drives the execution of the pseudo dojo tests."""
from __future__ import division, print_function

import sys

from argparse import ArgumentParser
from pprint import pprint

from pseudo_dojo import Dojo, TaskManager

__author__ = "Matteo Giantomassi"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"


def main():
    def str_examples():
        examples = """
          Usage Example:
          \n
          ppdojo_run.py Si.fhi -m10           => Cutoff converge test for Si.fhi using up to 10 CPUs.
          ppdojo_run.py Si.fhi -l1 -n10 -m20  => Deltafactor test using up to 20 CPUs, each task uses 10 MPI nodes.
        """
        return examples

    def show_examples_and_exit(error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        sys.exit(error_code)

    parser = ArgumentParser(epilog=str_examples())

    parser.add_argument('-l', '--max-level', type=int, default=0, 
                        help="Maximum DOJO level (default 0 i.e. ecut hints).")

    parser.add_argument('--loglevel', default="ERROR", type=str,
                         help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")

    parser.add_argument('pseudos', nargs='+', help='List of pseudopotential files.')

    try:
        options = parser.parse_args()
    except:
        shows_examples_and_exit(1)

    # loglevel is bound to the string value obtained from the command line argument. 
    # Convert to upper case to allow the user to specify --loglevel=DEBUG or --loglevel=debug
    import logging
    numeric_level = getattr(logging, options.loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % options.loglevel)
    logging.basicConfig(level=numeric_level)

    manager = TaskManager.from_user_config()
    pseudos = options.pseudos

    dojo = Dojo(pseudos[0], manager=manager, max_level=options.max_level)

    return dojo.start_training()

    #stats = []
    #for pseudo in pseudos:
    #    isok = dojo.add_pseudo(pseudo)
    #    stats.append(isok)

    #return stats.count(True)


if __name__ == "__main__":
    sys.exit(main())
