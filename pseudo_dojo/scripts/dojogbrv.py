#!/usr/bin/env python
"""Compute the deltafactor for a given pseudopotential."""
from __future__ import division, print_function, unicode_literals

import sys
import os
import argparse
import logging

from warnings import warn
from abipy import abilab
from pseudo_dojo.dojo.gbrv_outdb import GbrvOutdb
from pseudo_dojo.dojo.gbrv_compounds import GbrvCompoundsFactory

logger = logging.getLogger(__name__)


def gbrv_generate(options):
    """Generate the GBRV output databases."""

    # FIXME bugme with workdir
    for cls in GbrvOutdb.__subclasses__():
        outdb = cls.new_from_dojodir(options.dojo_dir)
        if os.path.exists(outdb.filepath):
            print("File %s already exists! " 
                  "New file won't be created. Remove it and try again"  % outdb.basename)
        else:
            outdb.json_write()
            print("Written new database %s" % outdb.basename)

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

    for path in options.database_list:
        outdb = GbrvOutdb.from_file(path)
        n = outdb.reset(status="failed")
        print("%s: %d has been resetted" % (outdb.basename, n))

    return 0


def gbrv_plot(options):
    """Plot data with matplotlib."""

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


def gbrv_run(options):
    options.manager = abilab.TaskManager.as_manager(options.manager)

    outdb = GbrvOutdb.from_file(options.database)

    jobs = outdb.find_jobs_torun(max_njobs=options.max_njobs)
    #jobs = outdb.find_jobs_torun(max_njobs=options.max_njobs, select_formulas=["NaCl"])
    num_jobs = len(jobs)

    if num_jobs == 0:
        print("Nothing to do, returning")
        return 0
    else:
        print("Will run %d works" % num_jobs)

    import tempfile
    workdir=tempfile.mkdtemp(dir=os.getcwd(), prefix=outdb.struct_type + "_")
    #workdir=tempfile.mkdtemp()

    flow = abilab.Flow(workdir=workdir, manager=options.manager)

    extra_abivars = {
        "mem_test": 0,
        "fband": 2,
        "nstep": 100,
        "paral_kgb": options.paral_kgb,
    }


    gbrv_factory = GbrvCompoundsFactory()
    for job in jobs:
        ecut = 30 if job.accuracy == "normal" else 45
        work = gbrv_factory.relax_and_eos_work(job.accuracy, job.pseudos, job.formula, outdb.struct_type, 
                                               ecut=ecut, pawecutdg=None, **extra_abivars)

        # Attach the database to the work to trigger the storage of the results.
        flow.register_work(work.set_outdb(outdb))

    print("Working in: ", flow.workdir)                     
    flow.build_and_pickle_dump()

    if options.dry_run:
        print("dry-run mode, will validate input files")
        isok, results = flow.abivalidate_inputs()              
        if not isok:                                           
            print(results)                                     
            return 1
    else:
        # Run the flow with the scheduler.
        #print("Fake Running")
        flow.use_smartio()
        flow.make_scheduler(rmflow=True).start()

    return 0


def main():
    def str_examples():
        return """\
usage example:
   dojogbrv generate directory      =>  Generate the json files needed GBRV computations.
                                        directory contains the pseudopotential table.
   dojogbrv update   directory      =>  Update all the json files in directory (check for 
                                        new pseudos or stale entries)
   dojogbrv reset dir/*.json        =>  Reset all failed entries in the json files
   dojogbrv run json_database       =>  Read data from json file, create flows and submit them
                                        with the scheduler.
"""

    def show_examples_and_exit(err_msg, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg: sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    # Parent parser for common options.
    copts_parser = argparse.ArgumentParser(add_help=False)
    copts_parser.add_argument('--loglevel', default="ERROR", type=str,
                               help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=str_examples(), formatter_class=argparse.RawDescriptionHelpFormatter)

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    # Subparser for the generate command.
    p_generate = subparsers.add_parser('generate', parents=[copts_parser], help="Denerate databases.")
    p_generate.add_argument('dojo_dir', help='Directory containing the pseudopotentials.')

    # Subparser for the update command.
    p_update = subparsers.add_parser('update', parents=[copts_parser], help="Update databases.")
    p_update.add_argument('dojo_dir', help='Directory containing the pseudopotentials.')

    # Subparser for the reset command.
    p_reset = subparsers.add_parser('reset', parents=[copts_parser], help="Reset failed entries in the database.")
    p_reset.add_argument('database_list', nargs="+", help='Database(s) with the output results.')

    # Subparser for plot command.
    p_plot = subparsers.add_parser('plot', parents=[copts_parser], help="Plot results in the databases.")
    p_plot.add_argument('database_list', nargs="+", help='Database(s) with the output results.')

    # Subparser for run command.
    p_run = subparsers.add_parser('run', parents=[copts_parser], help="Update databases.")
    p_run.add_argument('-m', '--manager', type=str, default=None,  help="Manager file")
    p_run.add_argument('--paral-kgb', type=int, default=0,  help="Paral_kgb input variable.")
    p_run.add_argument('-n', '--max-njobs', type=int, default=2, help="Maximum number of jobs (a.k.a. flows) that will be submitted")
    p_run.add_argument('-d', '--dry-run', default=False, action="store_true", help=("Dry run, build the flow and check validity "
                       "of input files without submitting"))
    p_run.add_argument('database', help='Database with the output results.')

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
    do_prof = False
    try:
        do_prof = sys.argv[1] == "prof"
        if do_prof: sys.argv.pop(1)
    except: 
        pass

    if do_prof:
        import pstats, cProfile
        cProfile.runctx("main()", globals(), locals(), "Profile.prof")
        s = pstats.Stats("Profile.prof")
        s.strip_dirs().sort_stats("time").print_stats()
    else:
        sys.exit(main())
