#!/usr/bin/env python
"""Script to run the GBRV tests for binary and ternary compunds."""
import sys
import os
import argparse
import hashlib

from monty.termcolor import cprint
from monty.functools import prof_main
from monty.io import FileLock
from pseudo_dojo.core.pseudos import DojoTable, OfficialDojoTable
from pseudo_dojo.refdata.gbrv.database import gbrv_database, species_from_formula
from pseudo_dojo.dojo.gbrv_outdb import GbrvOutdb
from pseudo_dojo.dojo.gbrv_compounds import GbrvCompoundsFactory, GbrvCompoundsFlow


def ecut_from_pseudos(pseudos, accuracy):
    """Compute ecut either from hints or from ppgen hints."""
    ecut, use_ppgen_hints = 0.0, False
    for p in pseudos:
        report = p.dojo_report
        if "hints" in report:
            ecut = max(ecut, report["hints"][accuracy]["ecut"])
        else:
            use_ppgen_hints = True
            ecut = max(ecut, report["ppgen_hints"][accuracy]["ecut"])

    assert ecut != 0.0
    if use_ppgen_hints:
        cprint("Hints are not available. Using ppgen_hints['high']", "yellow")

    return ecut


def gbrv_gendb(options):
    """Generate the GBRV output database from a djson file.."""
    # Build table from djson_path
    djson_path = os.path.abspath(options.djson_path)
    table = OfficialDojoTable.from_djson_file(djson_path)
    if options.verbose > 1: print(table)

    # Init database and dump it
    db = GbrvOutdb.new_from_table(table, djson_path)
    if os.path.exists(db.path):
        cprint("File %s already exists. New file won't be created. Remove it and try again" % db.path, "red")
        return 1

    db.json_write()
    cprint("Written new database %s" % os.path.relpath(db.path), "green")
    return 0


def gbrv_update(options):
    """Update the databases in dojo_dir."""
    raise NotImplementedError()
    filepath = os.path.join(options.dojo_dir, options.basename)
    outdb = GbrvOutdb.from_file(filepath)

    print("Checking:", outdb.basename)
    u = outdb.check_update()
    print("Update report:", u)

    return 0


def gbrv_rundb(options):
    """Build flow and run it."""
    dbpath = os.path.abspath(options.path)
    retcode = 0

    # Get list of jobs to execute.
    with FileLock(dbpath):
        outdb = GbrvOutdb.from_file(dbpath)
        jobs = outdb.find_jobs_torun(options.max_njobs)
        if not jobs:
            cprint("Nothing to do, returning 0", "yellow")
            return 0

    gbrv_factory = GbrvCompoundsFactory(xc=outdb["xc_name"])

    # Build workdir.
    s = "-".join(job.formula for job in jobs)
    m = hashlib.md5()
    m.update(s)
    workdir = os.path.join(os.getcwd(),
                          "GBRV_OUTDB_" + jobs[0].formula + "_" + jobs[-1].formula + "_" + m.hexdigest())
    #workdir = os.path.join(os.getcwd(), "GBRV_OUTDB_" + s)
    flow = GbrvCompoundsFlow(workdir=workdir)

    for job in jobs:
        #for accuracy in ("low", "normal", "high"):
        #for accuracy in ("high",):
        for accuracy in ("normal", "high"):
            ecut = max(p.hint_for_accuracy(accuracy).ecut for p in job.pseudos)
            pawecutdg = max(p.hint_for_accuracy(accuracy).pawecutdg for p in job.pseudos)
            if ecut <= 0.0: raise RuntimeError("Pseudos do not have hints")
            # Increase by 10 since many pseudos only have ppgen_hints
            #ecut += 10
            work = gbrv_factory.relax_and_eos_work(accuracy, job.pseudos, job.formula, job.struct_type,
                                                   ecut=ecut, pawecutdg=pawecutdg)

            # Attach the database to the work to trigger the storage of the results.
            flow.register_work(work.set_outdb(dbpath))

    print("Working in:", flow.workdir)
    flow.build_and_pickle_dump() #abivalidate=options.dry_run)
    if options.dry_run: return 0

    # Run the flow with the scheduler (enable smart_io)
    flow.use_smartio()
    retcode += flow.make_scheduler().start()

    return retcode


def gbrv_reset(options):
    """Reset entries in the databases."""
    status_list = []
    if "f" in options.status: status_list.append("failed")
    if "s" in options.status: status_list.append("scheduled")
    if not status_list:
        raise ValueError("Wrong value option %s" % options.status)
    cprint("Resetting all entries with status in: %s" % str(status_list), "yellow")

    outdb = GbrvOutdb.from_file(options.path)
    n = outdb.reset(status_list=status_list)
    print("%d entrie(s) have been resetted" % (n))

    return 0


def gbrv_plot(options):
    """Plot results in the database."""
    outdb = GbrvOutdb.from_file(options.path)
    frame = outdb.get_pdframe()
    #print(frame)
    #print(frame.describe())
    frame.print_summary()
    frame.plot_errors_for_elements()
    return 0

    # Use seaborn settings.
    import seaborn as sns
    sns.set(context=options.seaborn, style='darkgrid', palette='deep',
            font='sans-serif', font_scale=1, color_codes=False, rc=None)

    for struct_type in frame.struct_types():
        frame.plot_errors_for_structure(struct_type)
        frame.plot_hist(struct_type)

    return 0


def gbrv_notebook(options):
    """
    Generate ipython notebook to analyze the results in the database.
    """
    outdb = GbrvOutdb.from_file(options.path)
    return outdb.make_open_notebook()


def gbrv_runps(options):
    """
    Run GBRV compound tests given a list of pseudos.
    """
    # Build table and list of symbols
    pseudos = options.pseudos = DojoTable(options.pseudos)
    symbols = [p.symbol for p in pseudos]
    if options.verbose > 1: print(pseudos)

    # Consistency check
    assert len(set(symbols)) == len(symbols)
    xc_list = [p.xc for p in pseudos]
    xc = xc_list[0]
    if any(xc != other_xc for other_xc in xc_list):
        raise ValueError("Pseudos with different XC functional")

    gbrv_factory = GbrvCompoundsFactory(xc=xc)
    db = gbrv_factory.db

    entry = db.match_symbols(symbols)
    if entry is None:
        cprint("Cannot find entry for %s! Returning" % str(symbols), "red")
        return 1

    workdir = "GBRVCOMP_" + "_".join(p.basename for p in pseudos)
    print("Working in:", workdir)
    flow = GbrvCompoundsFlow(workdir=workdir)

    accuracy = "high"
    ecut = max(p.hint_for_accuracy(accuracy).ecut for p in pseudos)
    pawecutdg = max(p.hint_for_accuracy(accuracy).pawecutdg for p in pseudos)
    if ecut <= 0.0: raise RuntimeError("Pseudos do not have hints")
    #ecut = ecut_from_pseudos(pseudos, accuracy)
    print("Adding work for formula:", entry.symbol, ", structure:", entry.struct_type, ", ecut:", ecut)

    work = gbrv_factory.relax_and_eos_work(accuracy, pseudos, entry.symbol, entry.struct_type,
                                           ecut=ecut, pawecutdg=None)
    flow.register_work(work)

    flow.build_and_pickle_dump(abivalidate=options.dry_run)
    if options.dry_run: return 0

    # Run the flow with the scheduler (enable smart_io)
    flow.use_smartio()
    return flow.make_scheduler().start()


def gbrv_runform(options):
    """
    Run GBRV compound tests given a chemical formula.
    """
    # Extract chemical symbols from formula
    formula = options.formula
    symbols = set(species_from_formula(formula))

    # Init pseudo table and construct all possible combinations for the given formula.
    table = DojoTable.from_dir(top=options.pseudo_dir, exts=("psp8", "xml"), exclude_dirs="_*")
    pseudo_list = table.all_combinations_for_elements(symbols)

    #print("Removing relativistic pseudos from list")
    #pseudo_list = [plist for plist in pseudo_list if not any("_r" in p.basename for p in plist)]

    # This is hard-coded since we GBRV results are PBE-only.
    # There's a check between xc and pseudo.xc below.
    xc = "PBE"
    gbrv_factory = GbrvCompoundsFactory(xc=xc)
    db = gbrv_factory.db

    # Consistency check
    entry = db.match_symbols(symbols)
    if entry is None:
        cprint("Cannot find entry for %s! Returning" % str(symbols), "red")
        return 1

    workdir = "GBRVCOMP_" + formula
    print("Working in:", workdir)
    flow = GbrvCompoundsFlow(workdir=workdir)

    accuracy = "high"
    for pseudos in pseudo_list:
        if any(xc != p.xc for p in pseudos):
            raise ValueError("Pseudos with different XC functional")
        ecut = max(p.hint_for_accuracy(accuracy).ecut for p in pseudos)
        pawecutdg = max(p.hint_for_accuracy(accuracy).pawecutdg for p in pseudos)
        if ecut <= 0.0: raise RuntimeError("Pseudos do not have hints")
        #ecut = ecut_from_pseudos(pseudos, accuracy)
        print("Adding work for pseudos:", pseudos)
        print("    formula:", entry.symbol, ", structure:", entry.struct_type, ", ecut:", ecut)

        work = gbrv_factory.relax_and_eos_work(accuracy, pseudos, entry.symbol, entry.struct_type,
                                               ecut=ecut, pawecutdg=None)
        flow.register_work(work)

    if options.dry_run:
        flow.build_and_pickle_dump()
        return 0

    # Run the flow with the scheduler (enable smart_io)
    flow.use_smartio()
    return flow.make_scheduler().start()


def gbrv_find(options):
    """Print all formula containing symbols."""
    symbols = options.symbols
    print("Print all formula containing symbols: ", symbols)

    db = gbrv_database(xc="PBE")
    entries = db.entries_with_symbols(symbols)
    if not entries:
        cprint("Cannot find entries for %s! Returning" % str(symbols), "red")
        return 1

    print("Found %d entries" % len(entries))
    for i, entry in enumerate(entries):
        print("[%i]" % i, entry)

    return 0


def gbrv_info(options):
    """Print structure type and chemical formulas."""
    db = gbrv_database(xc="PBE")
    db.print_formulas()


@prof_main
def main():
    def str_examples():
        return """\
Usage example:
   dojogbrv.py info                         => Print all entries in the GBRV database.
   dojogbrv.py find Sr Si                   => Find entries containing these elements.
   dojogbrv.py runps Na/Na.psp8 F/F.psp8    => Run tests for a list of pseudos
   dojogbrv.py runform NaF -p pseudodir     => Run tests for NaF, take pseudos from pseudodir

   *** Under development. ***
   dojogbrv gendb table.djson               => Generate the json files needed for the GBRV computations.
                                               Requires path to djson file.
   dojogbrv rundb json_database             => Read data from json file, create flows and submit them
                                               with the scheduler.

   dojogbrv reset database --status=fs      => Reset all failed entries in the json database.
                                               f for failed, s for scheduled.
   dojogbrv update directory                => Update all the json files in directory (check for
                                               new pseudos or stale entries)
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
    parser = argparse.ArgumentParser(epilog=str_examples(), formatter_class=argparse.RawDescriptionHelpFormatter)

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    # Subparser for the gendb command.
    p_gendb = subparsers.add_parser('gendb', parents=[copts_parser], help=gbrv_gendb.__doc__)
    p_gendb.add_argument('djson_path', help='Path of the djson file.')

    # Subparser for the update command.
    p_update = subparsers.add_parser('update', parents=[copts_parser], help=gbrv_update.__doc__)
    p_update.add_argument('dojo_dir', help='Directory containing the pseudopotentials.')

    # Subparser for the reset command.
    p_reset = subparsers.add_parser('reset', parents=[copts_parser], help=gbrv_reset.__doc__)
    p_reset.add_argument('path', help='Database with the output results.')
    p_reset.add_argument("-s", '--status', type=str, default="f", help='f for failed, s for scheduled, `fs` for both')

    # Subparser for plot command.
    p_plot = subparsers.add_parser('plot', parents=[copts_parser], help=gbrv_plot.__doc__)
    p_plot.add_argument('path', help='path with the output results.')

    # Subparser for plot command.
    p_notebook = subparsers.add_parser('notebook', parents=[copts_parser], help=gbrv_notebook.__doc__)
    p_notebook.add_argument('path', help='path with the output results.')

    # Subparser for rundb command.
    p_rundb = subparsers.add_parser('rundb', parents=[copts_parser], help=gbrv_rundb.__doc__)

    #p_rundb.add_argument('--paral-kgb', type=int, default=0,  help="Paral_kgb input variable.")
    p_rundb.add_argument('-n', '--max-njobs', type=int, default=3,
                          help="Maximum number of jobs (a.k.a. works) that will be built and submitted")
    #def parse_formulas(s):
    #    return s.split(",") if s is not None else None
    #p_rundb.add_argument('-f', '--formulas', type=parse_formulas, default=None,
    #                    help="Optional list of formulas to be selected e.g. --formulas=LiF, NaCl")
    p_rundb.add_argument('path', help='path with the output results.')

    # Subparser for runps command.
    p_runps = subparsers.add_parser('runps', parents=[copts_parser], help=gbrv_runps.__doc__)
    p_runps.add_argument('pseudos', nargs="+", help="Pseudopotential files")

    p_runform = subparsers.add_parser('runform', parents=[copts_parser], help=gbrv_runform.__doc__)
    p_runform.add_argument('formula', help="Chemical formula.")
    p_runform.add_argument('-p', "--pseudo-dir", default=".", help="Directory with pseudos.")

    p_find = subparsers.add_parser('find', parents=[copts_parser], help=gbrv_find.__doc__)
    p_find.add_argument('symbols', nargs="+", help="Element symbols")

    p_info = subparsers.add_parser('info', parents=[copts_parser], help=gbrv_info.__doc__)

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
