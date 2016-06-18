#!/usr/bin/env python
"""Script to run the GBRV tests for binary and ternary compunds"""
from __future__ import division, print_function, unicode_literals

import sys
import os
import argparse
import logging

from monty.termcolor import cprint
from monty.functools import prof_main
from warnings import warn
from abipy import abilab
from pseudo_dojo.core.pseudos import DojoTable #dojopseudo_from_file,
from pseudo_dojo.refdata.gbrv import gbrv_database
from pseudo_dojo.dojo.gbrv_outdb import GbrvOutdb
from pseudo_dojo.dojo.gbrv_compounds import GbrvCompoundsFactory

logger = logging.getLogger(__name__)


def gbrv_dbgen(options):
    """Generate the GBRV output databases."""
    # FIXME bugme with workdir
    for cls in GbrvOutdb.__subclasses__():
        outdb = cls.new_from_dojodir(options.dojo_dir)
        if os.path.exists(outdb.filepath):
            cprint("File %s already exists! "
                  "New file won't be created. Remove it and try again"  % outdb.basename, "red")
        else:
            outdb.json_write()
            cprint("Written new database %s" % outdb.basename, "green")

    return 0


def gbrv_update(options):
    """Update the databases in dojo_dir."""
    for cls in GbrvOutdb.__subclasses__():
        filepath = os.path.join(options.dojo_dir, cls.basename)
        if not os.path.exists(filepath): continue

        outdb = cls.from_file(filepath)

        print("Checking:", outdb.basename)
        u = outdb.check_update()
        print("Update report:", u)

    return 0


def gbrv_reset(options):
    """Reset the failed entries in the list of databases specified by the user."""
    status_list = []
    if "f" in options.status: status_list.append("failed")
    if "s" in options.status: status_list.append("scheduled")
    if not status_list:
        raise ValueError("Wrong value option %s" % options.status)

    for path in options.database_list:
        outdb = GbrvOutdb.from_file(path)
        n = outdb.reset(status_list=status_list)
        print("%s: %d has been resetted" % (outdb.basename, n))
        outdb.json_write()

    return 0


def gbrv_plot(options):
    """Plot results in the databases."""

    for path in options.database_list:
        outdb = GbrvOutdb.from_file(path)
        frame = outdb.get_dataframe()
        print(frame)
        #import matplotlib.pyplot as plt
        #frame.plot(frame.index, ["normal_rel_err", "high_rel_err"])
        #ax.set_xticks(range(len(data.index)))
        #ax.set_xticklabels(data.index)
        #plt.show()
        #outdb.plot_errors(reference="ae", accuracy="normal")
        #for formula, records in outdb.values()
        #records = outdb["NaCl"]
        #for rec in records:
        #    rec.plot_eos()

    return 0


def gbrv_dbrun(options):
    """Build flow and run it."""
    outdb = GbrvOutdb.from_file(options.database)
    jobs = outdb.find_jobs_torun(max_njobs=options.max_njobs, select_formulas=options.formulas)
    num_jobs = len(jobs)
    if num_jobs == 0:
        cprint("Nothing to do, returning", "yellow")
        return 0
    else:
        print("Will run %d works" % num_jobs)

    gbrv_factory = GbrvCompoundsFactory(xc="PBE")

    import tempfile
    workdir=tempfile.mkdtemp(dir=os.getcwd(), prefix=outdb.struct_type + "_")
    #workdir=tempfile.mkdtemp()

    flow = abilab.Flow(workdir=workdir)

    extra_abivars = {
        "mem_test": 0,
        "fband": 2,
        "nstep": 100,
        "paral_kgb": options.paral_kgb,
    }

    for job in jobs:
        # FIXME this should be taken from the pseudos
        ecut = 30 if job.accuracy == "normal" else 45
        work = gbrv_factory.relax_and_eos_work(job.accuracy, job.pseudos, job.formula, outdb.struct_type,
                                               ecut=ecut, pawecutdg=None, **extra_abivars)

        # Attach the database to the work to trigger the storage of the results.
        flow.register_work(work.set_outdb(outdb.filepath))

    print("Working in: ", flow.workdir)
    flow.build_and_pickle_dump(abivalidate=options.dry_run)

    if not options.dry_run:
        # Run the flow with the scheduler (enable smart_io)
        flow.use_smartio()
        flow.make_scheduler().start()

    return 0


def gbrv_pprun(options):
    """
    Run GBRV compound tests for these pseudos.
    """
    pseudos = options.pseudos = DojoTable(options.pseudos)
    symbols = [p.symbol for p in pseudos]
    if options.verbose > 1: print(pseudos)

    xc_list = [p.xc for p in pseudos]
    xc = xc_list[0]
    if any(xc != other_xc for other_xc in xc_list):
        raise ValueError("Pseudos with different XC functional")

    gbrv_factory = GbrvCompoundsFactory(xc=xc)
    db = gbrv_factory.db
    #self.db = gbrv_database(xc)
    #for pseudo in pseudos:
    #    print(pseudo.symbol)
    #    for entry in db.entries_with_symbol(pseudo.symbol):
    #        print(entry)

    entries = db.entries_with_symbol_list(symbols)
    if not entries:
        cprint("Cannot find entries for %s! Returning" % str(symbols), "red")
        return 1
    #print(entries)

    workdir = "GBRV_COMPOUND_" + "_".join(p.basename for p in pseudos)
    print("Working in:", workdir)
    flow = abilab.Flow(workdir=workdir)
    #outdb = os.path.join(flow.workdir, "gbrv.json")

    # Get ecut from ppgen hints.
    ecut = 0.0
    for p in pseudos:
        # TODO: Should consider hints as well
        ecut = max(ecut, p.dojo_report["ppgen_hints"]["high"]["ecut"])
    assert ecut != 0.0

    for entry in entries:
        print("Adding work for formula:", entry.symbol, ", structure:", entry.struct_type, ", ecut:", ecut)
        work = gbrv_factory.relax_and_eos_work("normal", pseudos, entry.symbol, entry.struct_type,
                                               ecut=ecut, pawecutdg=None)
        #work.set_outdb(outdb)
        flow.register_work(work)

    flow.build_and_pickle_dump() #abivalidate=True)

    # Create output database.
    #import json
    #with open(outdb, "wt") as fh:
    #    print("Creating %s" % outdb)
    #    json.dump({}, fh)

    # Run the flow with the scheduler (enable smart_io)
    flow.use_smartio()
    flow.make_scheduler().start()


@prof_main
def main():
    def str_examples():
        return """\
usage example:
   dojogbrv dbgen directory      => Generate the json files needed for the GBRV computations.
                                    directory contains the pseudopotential table.
   dojogbrv update directory     =>  Update all the json files in directory (check for
                                     new pseudos or stale entries)
   dojogbrv reset dir/*.json     =>  Reset all failed entries in the json files
   dojogbrv dbrun json_database  =>  Read data from json file, create flows and submit them
                                     with the scheduler.
   dojogbrv pprun pseudos        =>  create flows and submit them with the scheduler.
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

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=str_examples(), formatter_class=argparse.RawDescriptionHelpFormatter)

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    # Subparser for the dbgen command.
    p_dbgen = subparsers.add_parser('dbgen', parents=[copts_parser], help=gbrv_dbgen.__doc__)
    p_dbgen.add_argument('dojo_dir', help='Directory containing the pseudopotentials.')

    # Subparser for the update command.
    p_update = subparsers.add_parser('update', parents=[copts_parser], help=gbrv_update.__doc__)
    p_update.add_argument('dojo_dir', help='Directory containing the pseudopotentials.')

    # Subparser for the reset command.
    p_reset = subparsers.add_parser('reset', parents=[copts_parser], help=gbrv_reset.__doc__)
    p_reset.add_argument("-s", '--status', type=str, default="f", help='f for failed, s for scheduled, `fs` for both')
    p_reset.add_argument('database_list', nargs="+", help='Database(s) with the output results.')

    # Subparser for plot command.
    p_plot = subparsers.add_parser('plot', parents=[copts_parser], help=gbrv_plot.__doc__)
    p_plot.add_argument('database_list', nargs="+", help='Database(s) with the output results.')

    # Subparser for dbrun command.
    p_dbrun = subparsers.add_parser('dbrun', parents=[copts_parser], help=gbrv_dbrun.__doc__)

    p_dbrun.add_argument('--paral-kgb', type=int, default=0,  help="Paral_kgb input variable.")
    p_dbrun.add_argument('-n', '--max-njobs', type=int, default=2,
                          help="Maximum number of jobs (a.k.a. flows) that will be submitted")
    p_dbrun.add_argument('-d', '--dry-run', default=False, action="store_true",
                         help=("Dry run, build the flow and check validity of input files without submitting"))
    def parse_formulas(s):
        return s.split(",") if s is not None else None
    p_dbrun.add_argument('-f', '--formulas', type=parse_formulas, default=None,
                        help="Optional list of formulas to be selected e.g. --formulas=LiF, NaCl")
    p_dbrun.add_argument('database', help='Database with the output results.')

    # Subparser for pprun command.
    p_pprun = subparsers.add_parser('pprun', parents=[copts_parser], help=gbrv_pprun.__doc__)
    p_pprun.add_argument('pseudos', nargs="+", help="Pseudopotential files")

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

    # Dispatch.
    return globals()["gbrv_" + options.command](options)


if __name__ == "__main__":
    sys.exit(main())
