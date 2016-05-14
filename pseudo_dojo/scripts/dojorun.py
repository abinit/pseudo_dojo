#!/usr/bin/env python
"""Compute the deltafactor for a given pseudopotential."""
from __future__ import division, print_function, unicode_literals

import sys
import os
import argparse
import copy
import numpy as np
import logging
import abipy.abilab as abilab

from warnings import warn
from pseudo_dojo.dojo.works import DeltaFactory, GbrvFactory, DFPTPhononFactory, EbandsFactory
from pymatgen.io.abinit.pseudos import Pseudo
from pymatgen.core.periodic_table import PeriodicTable


logger = logging.getLogger(__name__)


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
    try:
        ppgen_ecut = int(report["ppgen_hints"]["high"]["ecut"])
        ecut_list = copy.copy(report["ecuts"])

    except KeyError:
        print('New pseudo without report from the generator, the convergence study is started from 16H')
        report["ppgen_hints"] = {}
        report["ppgen_hints"]["high"] = {} 
        report["ppgen_hints"]["high"]["ecut"] = 16.0
        report["ecuts"] = [16.0, 20.0, 24.0]
        pseudo.write_dojo_report(report)
        ppgen_ecut = int(report["ppgen_hints"]["high"]["ecut"])
        ecut_list = copy.copy(report["ecuts"])

    try:
        ecut_hint = int(report["hints"]["normal"]["ecut"])
    except KeyError:
        try:
            ecut_hint = int(report["ppgen_hints"]["normal"]["ecut"])
        except KeyError:
            ecut_hint = ppgen_ecut

    #if 'extend' in options:
    #    next_ecut = max(ecut_list) + 2
    #    ecut_list.append(next_ecut)

    #if 'new-ecut' in options:
    #    ecut_list.append(options['new-ecut'])

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
        #FIXME
        #factory = DeltaFactory(xc=pseudo.xc)
        if os.path.isfile('LDA'):
            factory = DeltaFactory(xc='LDA')
        else:
            factory = DeltaFactory()

        for ecut in ecut_list:
            str_ecut = '%.1f' % ecut
            if "deltafactor" in report and str_ecut in report["deltafactor"]:
                print("[deltafactor]: ignoring ecut=", str_ecut, "because it's already in the DOJO_REPORT")
                continue

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
                str_ecut = '%.1f' % ecut
                if dojo_trial in report and str_ecut in report[dojo_trial]:
                    print("[gbrv]: ignoring ecut=", str_ecut, "because it's already in the DOJO_REPORT")
                    continue

                pawecutdg = 2 * ecut if pseudo.ispaw else None
                # FIXME: we use ntime=3, because structure relaxations go bananas after the third step.
                work = gbrv_factory.relax_and_eos_work(pseudo, struct_type, ecut=ecut, ntime=5, pawecutdg=pawecutdg, **extra_abivars)
                flow.register_work(work, workdir="GBRV_" + struct_type + str(ecut))

    # PHONON test
    if "phonon" in options.trials:
        phonon_factory = DFPTPhononFactory()

        for ecut in ecut_list:
            str_ecut = '%.1f' % ecut
            if "phonon" in report and str_ecut in report["phonon"]:
                print("[phonon]: ignoring ecut=", str_ecut, "because it's already in the DOJO_REPORT")
                continue

            kppa = 1000
            pawecutdg = 2 * ecut if pseudo.ispaw else None
            work = phonon_factory.work_for_pseudo(pseudo, accuracy="high", kppa=kppa, ecut=ecut, pawecutdg=pawecutdg,
                                                  tolwfr=1.e-20, smearing="fermi_dirac:0.0005", qpt=[0,0,0], mem_test=0)
            if work is not None:
                flow.register_work(work, workdir='GammaPhononsAt'+str(ecut))
            else:
                warn('cannot create GammaPhononsAt' + str(ecut) + ' work, factory returned None')

    # PHONON WihtOut Asr test
    if "phwoa" in options.trials:
        phonon_factory = DFPTPhononFactory()

        for ecut in [ecut_list[0], ecut_list[-1]]:
            str_ecut = '%.1f' % ecut
            if "phwoa" in report and str_ecut in report["phwoa"]:
                print("[phwoa]: ignoring ecut=", str_ecut, "because it's already in the DOJO_REPORT")
                continue

            kppa = 1000
            pawecutdg = 2 * ecut if pseudo.ispaw else None
            work = phonon_factory.work_for_pseudo(pseudo, accuracy="high", kppa=kppa, ecut=ecut, pawecutdg=pawecutdg,
                                                  tolwfr=1.e-20, smearing="fermi_dirac:0.0005", qpt=[0,0,0], rfasr=0)
            if work is not None:
                flow.register_work(work, workdir='GammaPhononsAt'+str(ecut)+'WOA')
            else:
                warn('cannot create GammaPhononsAt' + str(ecut) + 'WOA work, factory returned None')

    # EBANDS test
    if "ebands" in options.trials:
        ebands_factory = EbandsFactory()
        ecut = ecut_hint    
        str_ecut = '%.1f' % ecut

        if "ebands" in report and str_ecut in report["ebands"]:
            print("[ebands]: ignoring ecut=", str_ecut, "because it's already in the DOJO_REPORT")
        else:
            kppa = 3000
            pawecutdg = 2 * ecut if pseudo.ispaw else None
            work = ebands_factory.work_for_pseudo(pseudo, accuracy="high", kppa=kppa, ecut=ecut, pawecutdg=pawecutdg,
                                                  bands_factor=15, smearing="fermi_dirac:0.0005", qpt=[0,0,0], mem_test=0)
            if work is not None:
                flow.register_work(work, workdir='EbandsAt' + str(ecut))
            else:
                warn('cannot create EbandsAt' + str(ecut) + ' work, factory returned None')

    if len(flow) > 0:
        return flow.allocate()
    else:
        # Empty flow since all trials have been already performed.
        return None


def main():
    def str_examples():
        examples = """
                   Usage Example:\n
                   ppdojo_run.py Si.psp8  => Build pseudo_dojo flow for Si.fhi\n
                   """
        return examples

    def show_examples_and_exit(error_code=1):
        """Display the usage of the script."""
        print(str_examples())
        sys.exit(error_code)

    parser = argparse.ArgumentParser(epilog=str_examples())

    parser.add_argument('-m', '--manager', type=str, default=None,  help="Manager file")
    parser.add_argument('-d', '--dry-run', default=False, action="store_true", help="Dry run, build the flow without submitting it")
    parser.add_argument('--paral-kgb', type=int, default=0,  help="Paral_kgb input variable.")
    parser.add_argument('-p', '--plot', default=False, action="store_true", help="Plot convergence when the flow is done")
    parser.add_argument('-n', '--new-ecut', type=int, default=None, action="store", help="Extend the ecut grid with the new-ecut point")

    def parse_trials(s):
        if s == "all": return ["df", "gbrv", "phonon", "phowa"]
        return s.split(",")

    parser.add_argument('--trials', default="all",  type=parse_trials, help=("List of tests e.g --trials=df,gbrv,phonon,phwoa\n"
                        "  df:     test delta factor against all electron refference\n"
                        "  gbrv:   test fcc and bcc lattice parameters agains AE refference\n"
                        "  phonon: test phonon mode at gamma convergence\n"
                        "  phwoa:  test violation of the acoustic sum rule (without enforcing it) at the min and max ecut\n"))

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
            flow.make_scheduler().start()

    else:
        # Gather all pseudos starting from the current working directory and run the flows iteratively.
        table = PeriodicTable()
        all_symbols = set(element.symbol for element in table.all_elements)
        #all_symbols = ["H"]
        #print(os.listdir(options.path))

        #print("here", os.path.basename(os.path.dirname(options.path)))
        #print("here", options.path)
        if os.path.basename(os.path.dirname(options.path)) in all_symbols:
            #print("here")
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

        print(pseudos)

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

            new_report = pseudo.read_dojo_report()
            new_report.plot_deltafactor_convergence()
            new_report.plot_gbrv_convergence()
            new_report.plot_phonon_convergence()

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
