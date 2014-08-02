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


def str_examples():
    examples = """
      Usage Example:
      \n
      ppdojo_run.py Si.fhi -m10           => Cutoff converge test for Si.fhi using up to 10 CPUs.
      ppdojo_run.py Si.fhi -l1 -n10 -m20  => Deltafactor test using up to 20 CPUs, each task uses 10 MPI nodes.
    """
    return examples

def show_examples_and_exit(err_msg=None, error_code=1):
    "Display the usage of the script."
    sys.stderr.write(str_examples())
    if err_msg: 
        sys.stderr.write("Fatal Error\n" + err_msg + "\n")
    sys.exit(error_code)


def main():
    parser = ArgumentParser(epilog=str_examples())

    parser.add_argument('-m', '--max_ncpus', type=int, default=1,
                        help="Maximum number of CPUs that can be used by the DOJO.")

    parser.add_argument('-n', '--mpi-ncpus', type=int, default=1,
                        help="Number of MPI processes per task).")

    parser.add_argument('-l', '--max-level', type=int, default=0, 
                        help="Maximum DOJO level (default 0 i.e. ecut hints).")

    parser.add_argument('-a', '--accuracy', type=str, default="normal", 
                        help="Accuracy of the calculation (low-normal-high). Default: normal.")


    parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
                         help='Verbose, can be supplied multiple times to increase verbosity')

    parser.add_argument('--loglevel', default="ERROR", type=str,
                         help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")

    parser.add_argument('pseudos', nargs='+', help='List of pseudopotential files.')

    options = parser.parse_args()

    # loglevel is bound to the string value obtained from the command line argument. 
    # Convert to upper case to allow the user to specify --loglevel=DEBUG or --loglevel=debug
    import logging
    numeric_level = getattr(logging, options.loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % options.loglevel)
    logging.basicConfig(level=numeric_level)

    pseudos = options.pseudos
    max_ncpus = options.max_ncpus
    mpi_ncpus = options.mpi_ncpus

    if mpi_ncpus > max_ncpus:
        err_msg = "mpi_cpus %(mpi_ncpus)s cannot be greater than max_ncpus %(max_ncpus)s" % locals()
        show_examples_and_exit(err_msg=err_msg)

    manager = TaskManager.from_user_config()
    print(manager)

    dojo = Dojo(manager=manager,
                max_ncpus=max_ncpus, 
                max_level=options.max_level, 
                verbose=options.verbose)

    stats = []
    for pseudo in pseudos:
        isok = dojo.challenge_pseudo(pseudo, accuracy=options.accuracy)
        stats.append(isok)

    return stats.count(True)


if __name__ == "__main__":
    sys.exit(main())
