#!/usr/bin/env python
"""This scripts drives the execution of the pseudo dojo tests."""
from __future__ import division, print_function

import sys

from argparse import ArgumentParser
from pprint import pprint

from pseudo_dojo.core import RunMode, Dojo

__author__ = "Matteo Giantomassi"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"

################################################################################

def str_examples():
    examples = """
      Example usage:
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

    parser.add_argument('pseudos', nargs='+', help='List of pseudopotential files.')

    options = parser.parse_args()

    pseudos = options.pseudos
    max_ncpus = options.max_ncpus
    mpi_ncpus = options.mpi_ncpus

    if mpi_ncpus > max_ncpus:
        err_msg = "mpi_cpus %(mpi_ncpus)s cannot be greater than max_ncpus %(max_ncpus)s" % locals()
        show_examples_and_exit(err_msg=err_msg)

    runmode = RunMode.mpi_parallel(mpi_ncpus=mpi_ncpus)
    #runmode = RunMode.load_user_configuration()
    pprint(runmode)

    dojo = Dojo(runmode=runmode, 
                max_ncpus=max_ncpus, 
                max_level=options.max_level, 
                verbose=options.verbose,
                )

    for pseudo in pseudos:
        dojo.challenge_pseudo(pseudo, accuracy=options.accuracy)

    return 0

################################################################################

if __name__ == "__main__":
    sys.exit(main())
