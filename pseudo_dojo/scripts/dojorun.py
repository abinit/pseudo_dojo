#!/usr/bin/env python
"""Compute the deltafactor for a given pseudopotential."""
from __future__ import division, print_function, unicode_literals

import sys
import os
import argparse
import numpy as np
import abipy.abilab as abilab

from warnings import warn
from pseudo_dojo.dojo.works import DeltaFactory, GbrvFactory
from pymatgen.io.abinitio.pseudos import Pseudo
from pymatgen.core.periodic_table import PeriodicTable


class RedirectStdStreams(object):
    """See http://stackoverflow.com/questions/6796492/temporarily-redirect-stdout-stderr"""
    def __init__(self, stdout=None, stderr=None):
        self._stdout = stdout or sys.stdout
        self._stderr = stderr or sys.stderr

    def __enter__(self):
        self.old_stdout, self.old_stderr = sys.stdout, sys.stderr
        self.old_stdout.flush(); self.old_stderr.flush()
        sys.stdout, sys.stderr = self._stdout, self._stderr

    def __exit__(self, exc_type, exc_value, traceback):
        self._stdout.flush(); self._stderr.flush()
        sys.stdout = self.old_stdout
        sys.stderr = self.old_stderr


def build_flow(pseudo, options):
    """Build the flow, returns None if no calculation must be performed.""" 
    pseudo = Pseudo.as_pseudo(pseudo)

    workdir = pseudo.basename + "_DOJO"
    if os.path.exists(workdir): 
        warn("Directory %s already exists" % workdir)
        return None
        #raise ValueError("%s exists" % workdir)

    flow = abilab.Flow(workdir=workdir, manager=options.manager)

    extra_abivars = {
            "mem_test": 0,
            "fband": 2,
            "nstep": 100,
            "paral_kgb": options.paral_kgb
            #"nsym": 1,
    }

    report = pseudo.read_dojo_report()
    #print(report)
    #hints = report["hints"]

    # Build ecut mesh.
    ppgen_ecut = int(report["ppgen_hints"]["high"]["ecut"])

    ecut_list = report["ecuts"]

    add_ecuts = False
    if add_ecuts:
        #dense_right = np.linspace(ppgen_ecut, ppgen_ecut + 10, num=6)
        #dense_left = np.linspace(ppgen_ecut-8, ppgen_ecut, num=4, endpoint=False)
        #coarse_high = np.linspace(ppgen_ecut + 15, ppgen_ecut + 40, num=4)

        dense_right = np.arange(ppgen_ecut, ppgen_ecut + 6*2, step=2)
        dense_left = np.arange(max(ppgen_ecut-6, 2), ppgen_ecut, step=2)
        coarse_high = np.arange(ppgen_ecut + 15, ppgen_ecut + 35, step=5)

        ecut_list = list(dense_left) + list(dense_right) + list(coarse_high)

    # Computation of the deltafactor.
    if "df" in options.trials:
        factory = DeltaFactory()
        for ecut in ecut_list:
            if "deltafactor" in report and ecut in report["deltafactor"].keys(): continue
            pawecutdg = 2 * ecut if pseudo.ispaw else None
            # Build and register the workflow.
            work = factory.work_for_pseudo(pseudo, kppa=6750, ecut=ecut, pawecutdg=pawecutdg, **extra_abivars)
            flow.register_work(work, workdir='WDF' + str(ecut))

    # GBRV tests.
    if "gbrv" in options.trials:
        gbrv_factory = GbrvFactory()
        gbrv_structs = ("fcc", "bcc")
        for struct_type in gbrv_structs:
            dojo_trial = "gbrv_" + struct_type
            for ecut in ecut_list:
                if dojo_trial in report and ecut in report[dojo_trial].keys(): continue
                pawecutdg = 2 * ecut if pseudo.ispaw else None
                work = gbrv_factory.relax_and_eos_work(pseudo, struct_type, ecut=ecut, pawecutdg=pawecutdg, **extra_abivars)
                flow.register_work(work, workdir="GBRV_" + struct_type + str(ecut))

    if len(flow) > 0:
        return flow.allocate()
    else:
        # Empty flow since all trials have been already performed.
        return None


def main():
    def str_examples():
        examples = """
Usage Example:\n
    ppdojo_run.py Si.psp8  => Build pseudo_dojo flow for Si.fhi
\n"""
        return examples

    def show_examples_and_exit(error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        sys.exit(error_code)

    parser = argparse.ArgumentParser(epilog=str_examples())

    parser.add_argument('-m', '--manager', type=str, default=None,  help="Manager file")
    parser.add_argument('-d', '--dry-run', type=bool, default=False,  help="Dry run, build the flow without submitting it")
    parser.add_argument('--paral-kgb', type=int, default=1,  help="Paral_kgb input variable.")

    def parse_trials(s):
        if s == "all": return ["df", "gbrv"]
        return s.split(",")

    parser.add_argument('--trials', default="all",  type=parse_trials, help="List of tests e.g --trials=df,gbrv")

    parser.add_argument('--loglevel', default="ERROR", type=str,
                        help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")

    parser.add_argument('path', help='pseudopotential file.')

    # Create the parsers for the sub-commands
    #subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")
    # Subparser for single command.
    #p_build = subparsers.add_parser('build', help="Build dojo.")

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

    options.manager = abilab.TaskManager.from_user_config() if options.manager is None else \
                      abilab.TaskManager.from_file(options.manager)

    if os.path.isfile(options.path):
        # Operate on a single pseudo.
        flow = build_flow(options.path, options)
        if flow is None: 
            warn("DOJO_REPORT is already computed for pseudo %s." % options.path)
            return 0
        if options.dry_run:
            flow.build_and_pickle_dump()
        else:
            # Run the flow with the scheduler.
            #print("nlaunch: %d" % flow.rapidfire())
            flow.make_scheduler().start()

    else:
        # Gather all pseudos starting from the current working directory and run the flows iteratively.
        table = PeriodicTable()
        all_symbols = set(element.symbol for element in table.all_elements)
        #all_symbols = ["H"]
        #print(os.listdir(options.path))

        print("here" , os.path.basename(os.path.dirname(options.path)))
        print("here" , options.path)
        if os.path.basename(os.path.dirname(options.path)) in all_symbols:
            print("here")
            dirs = [options.path]
        else:
            dirs = [os.path.join(options.path, d) for d in os.listdir(options.path) if d in all_symbols]
        print(dirs)

        pseudos = []
        for d in dirs:
            #print(d)
            pseudos.extend(os.path.join(d, p) for p in os.listdir(d) if p.endswith(".psp8"))

        if not pseudos:
            warn("Empty list of pseudos")
            return 0

        nflows, nlaunch = 0, 0
        #exc_filename = "allscheds_exceptions.log"
        #if os.path.exists(exc_filename):
        #    raise RuntimeError("File %s already exists, remove it before running the script" % exc_filename)
        #exc_log = open(exc_filename, "w")
        exc_log = sys.stderr

        for pseudo in pseudos:
            pseudo = Pseudo.as_pseudo(pseudo)
            report = pseudo.dojo_report
            if "version" not in report: continue

            flow = build_flow(pseudo, options)
            if flow is None: 
                warn("DOJO_REPORT is already computed for pseudo %s." % pseudo.basename)
                continue

            #if os.path.exists(flow.workdir) or nflows >= 2: continue
            nflows += 1

            try:
                flow.make_scheduler().start()
            except Exception as exc:
                # Log exception and proceed with the next pseudo.
                exc_log.write(str(exc))

            #with open(pseudo.basename + "sched.stdout", "w") as sched_stdout, \
            #     open(pseudo.basename + "sched.stderr", "w") as sched_stderr: 
            #    with RedirectStdStreams(stdout=sched_stdout, stderr=sched_stderr):
            #        try:
            #            flow.make_scheduler().start()
            #        except Exception as exc:
            #            # Log exception and proceed with the next pseudo.
            #            exc_log.write(str(exc))

        #exc_log.close()
        #print("nlaunch: %d" % nlaunch)
        #print("nflows: %d" % nflows)

    return 0


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
