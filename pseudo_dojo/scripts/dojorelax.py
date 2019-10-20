#!/usr/bin/env python
"""Script to perform structural relaxation."""
import sys
import os
import argparse
import logging

from monty.termcolor import cprint
from monty.functools import prof_main
from abipy import abilab
from pseudo_dojo.core.pseudos import dojopseudo_from_file

logger = logging.getLogger(__name__)


def ecut_from_pseudo(pseudo):
    """Compute ecut either from hints or from ppgen hints."""
    ecut, use_ppgen_hints = 0.0, False
    report = pseudo.dojo_report
    if "hints" in report:
        ecut = max(ecut, report["hints"]["high"]["ecut"])
    else:
        use_ppgen_hints = True
        ecut = max(ecut, report["ppgen_hints"]["high"]["ecut"])

    assert ecut != 0.0
    if use_ppgen_hints:
        cprint("Hints are not available. Using ppgen_hints['high']", "yellow")
        #ecut += 10

    return ecut


@prof_main
def main():
    def str_examples():
        return """\
Usage example:
   dojogbrv.py info                         => Print all entries in the GBRV database.
"""

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg: sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    # Parent parser for common options.
    copts_parser = argparse.ArgumentParser(add_help=False)
    copts_parser.add_argument('--loglevel', default="ERROR", type=str,
                              help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")
    copts_parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
                              help='Verbose, can be supplied multiple times to increase verbosity')
    copts_parser.add_argument('-d', '--dry-run', default=False, action="store_true",
                              help=("Dry run, build the flow and check validity of input files without submitting"))

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=str_examples(), parents=[copts_parser],
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('pseudo', help='Pseudopotential file.')

    # Create the parsers for the sub-commands
    #subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    # Subparser for the gendb command.
    #p_gendb = subparsers.add_parser('gendb', parents=[copts_parser], help=gbrv_gendb.__doc__)
    #p_gendb.add_argument('djson_path', help='Path of the djson file.')

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

    pseudo = options.pseudo = dojopseudo_from_file(options.pseudo)
    ecut = ecut_from_pseudo(pseudo)

    workdir = pseudo.basename + "_DOJO_RELAX"
    if os.path.exists(workdir):
        cprint("Directory %s already exists" % workdir, "red")
        return 1

    flow = abilab.Flow(workdir=workdir)
    pawecutdg = 2 * ecut if pseudo.ispaw else None

    jth_afcc = dict(
        La=5.272,
        Ce=4.769,
        Pr=4.608,
        Nd=4.531,
        Pm=4.499,
        Sm=4.516,
        Eu=4.630,
        Gd=4.723,
        Tb=4.835,
        Dy=4.918,
        Ho=5.012,
        Er=5.073,
        Tm=5.120,
        Yb=5.163,
        Lu=5.589,
    )

    a_guess = jth_afcc[pseudo.symbol]

    from pseudo_dojo.dojo.works import RelaxWithGbrvParamsWork
    ecut_list = [ecut, ecut+10, ecut+20]

    work = RelaxWithGbrvParamsWork(
                 a_guess, "fcc", pseudo, ecut_list=ecut_list, pawecutdg=None,
                 ngkpt=(8, 8, 8))
                 #ngkpt=(12, 12, 12))
                 #spin_mode="unpolarized", include_soc=False, tolvrs=1.e-10, smearing="fermi_dirac:0.001 Ha",
                 #ecutsm=0.05, chksymbreak=0)
    flow.register_work(work)

    #flow.build_and_pickle_dump(abivalidate=True)
    #return 0
    return flow.make_scheduler().start()

    # Dispatch.
    #return globals()["gbrv_" + options.command](options)


if __name__ == "__main__":
    sys.exit(main())
