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

from monty.termcolor import cprint
from monty.functools import prof_main
from pymatgen.io.abinit.pseudos import Pseudo
from pymatgen.core.periodic_table import PeriodicTable
from pseudo_dojo.core.pseudos import dojopseudo_from_file
from pseudo_dojo.dojo.works import DeltaFactory, GbrvFactory, DFPTPhononFactory, EbandsFactory


logger = logging.getLogger(__name__)


def build_flow(pseudo, options):
    """Build the flow, returns None if no calculation must be performed."""
    print(pseudo)
    if not pseudo.has_dojo_report:
        raise ValueError("Cannot find dojo_report")

    if options.soc and not pseudo.supports_soc:
        raise TypeError("SOC is on but pseudo does not support spin-orbit coupling")

    if not options.soc and pseudo.supports_soc and pseudo.path.endswith("psp8"):
        cprint("[STRANGE]: Your psp8 pseudo supports SOC but options.soc is off", "magenta")

    report = pseudo.dojo_report
    if options.verbose > 1: print(report)

    workdir = pseudo.basename + "_DOJO"
    if os.path.exists(workdir):
        cprint("Directory %s already exists" % workdir, "red")
        return None

    flow = abilab.Flow(workdir=workdir, manager=options.manager)

    extra_abivars = {
            "mem_test": 0,
            #"nstep": 100,
            "paral_kgb": options.paral_kgb
    }
    #flow.walknset_vars(extra_abivars)

    # Build ecut mesh.
    try:
        ppgen_ecut = int(report["ppgen_hints"]["high"]["ecut"])
        ecut_list = copy.copy(report["ecuts"])

    except KeyError:
        cprint('New pseudo without report from the generator, the convergence study is started from 16H', "yellow")
        #raise NotImplementedError()
        # TODO
        #report = DojoReport.from_pseudo(pseudo)
        report["ppgen_hints"] = {}
        report["ppgen_hints"]["high"] = {}
        report["ppgen_hints"]["high"]["ecut"] = 16.0
        report["ecuts"] = [16.0, 20.0, 24.0]
        report.json_write(pseudo.djrepo_path)
        #pseudo.write_dojo_report(report)
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
        factory = DeltaFactory(xc=pseudo.xc)
        dojo_trial = "deltafactor" is not options.soc else "deltafactor_soc"

        for ecut in ecut_list:
            if report.has_trial(dojo_trial, ecut=ecut):
                cprint("[%s]: ignoring ecut=%s because it's already in the DOJO_REPORT" % (dojo_trial, ecut), "magenta")
                continue

            # Build and register the work.
            pawecutdg = 2 * ecut if pseudo.ispaw else None
            work = factory.work_for_pseudo(pseudo, kppa=6750, ecut=ecut, pawecutdg=pawecutdg,
                                           include_soc=options.soc, **extra_abivars)
            flow.register_work(work, workdir='WDF' + str(ecut))

    # GBRV tests.
    if "gbrv" in options.trials:
        gbrv_factory = GbrvFactory(xc=pseudo.xc)
        gbrv_structs = ("fcc", "bcc")

        for struct_type in gbrv_structs:
            dojo_trial = "gbrv_" + struct_type
            for ecut in ecut_list:
                if report.has_trial(dojo_trial, ecut=ecut):
                    cprint("[gbrv]: ignoring ecut=%s because it's already in the DOJO_REPORT" % ecut, "magenta")
                    continue

                # Build and register the work.
                pawecutdg = 2 * ecut if pseudo.ispaw else None
                work = gbrv_factory.relax_and_eos_work(pseudo, struct_type, ecut=ecut, pawecutdg=pawecutdg,
                                                       include_soc=options.soc, ntime=50, **extra_abivars)
                flow.register_work(work, workdir="GBRV_" + struct_type + str(ecut))

    # PHONON test
    if "phonon" in options.trials:
        assert not options.soc
        phonon_factory = DFPTPhononFactory(xc=pseudo.xc)
        dojo_trial = "phonon"

        for ecut in ecut_list:
            if report.has_trial(dojo_trial, ecut=ecut):
                cprint("[phonon]: ignoring ecut=%s because it's already in the DOJO_REPORT" % ecut, "magenta")
                continue

            # Build and register the work.
            pawecutdg = 2 * ecut if pseudo.ispaw else None
            work = phonon_factory.work_for_pseudo(pseudo, kppa=1000, ecut=ecut, pawecutdg=pawecutdg,
                                                  tolwfr=1.e-20, smearing="fermi_dirac:0.0005", qpt=[0,0,0], mem_test=0)
            if work is not None:
                flow.register_work(work, workdir='GammaPhononsAt' + str(ecut))
            else:
                cprint('Cannot create GammaPhononsAt%s work, factory returned None' % str(ecut), "magenta")

    # PHONON WihtOut Asr test
    if "phwoa" in options.trials:
        assert not options.soc
        phonon_factory = DFPTPhononFactory(xc=pseudo.xc)

        for ecut in [ecut_list[0], ecut_list[-1]]:
            str_ecut = '%.1f' % ecut
            if "phwoa" in report and str_ecut in report["phwoa"]:
                cprint("[phwoa]: ignoring ecut=%s because it's already in the DOJO_REPORT" % str_ecut, "magenta")
                continue

            pawecutdg = 2 * ecut if pseudo.ispaw else None
            work = phonon_factory.work_for_pseudo(pseudo, kppa=1000, ecut=ecut, pawecutdg=pawecutdg,
                                                  tolwfr=1.e-20, smearing="fermi_dirac:0.0005", qpt=[0,0,0], rfasr=0)
            if work is not None:
                flow.register_work(work, workdir='GammaPhononsAt'+str(ecut)+'WOA')
            else:
                cprint('Cannot create GammaPhononsAt %s WOA work, factory returned None' % str(ecut), "red")

    # GHOSTS test
    if "ghosts" in options.trials:
        assert not options.soc
        dojo_trial = "ghosts" is not options.soc else "ghosts_soc"
        ebands_factory = EbandsFactory(pseudo.xc)
        ecut = int(report["ppgen_hints"]["high"]["ecut"])
        pawecutdg = None if not pseudo.ispaw else int(report["ppgen_hints"]["high"]["pawecutdg"])

        #str_ecut = '%.1f' % ecut
        #print(report["ghosts"].pop(ecut, None))

        if report.has_trial(dojo_trial, ecut=ecut):
            cprint("[%s]: ignoring ecut=%s because it's already in the DOJO_REPORT" % (dojo_trial, ecut), "magenta")
        else:
            # Build and register the work.
            work = ebands_factory.work_for_pseudo(pseudo, kppa=2000, maxene=250,
                                                  ecut=ecut, pawecutdg=pawecutdg,
                                                  **extra_abivars)
            if work is not None:
                flow.register_work(work, workdir='EbandsAt' + str(ecut))
            else:
                cprint('Cannot create EbandsAt%s work, factory returned None' % str(ecut), "magenta")

    # phonons at gamma test.
    if "phgamma" in options.trials:
        assert not options.soc
        phg_factory = GammaPhononFactory(pseudo.xc)
        dojo_trial = "phgamma" if not options.soc else "phgamma_soc"

        for ecut in ecut_list:
            if report.has_trial(dojo_trial, ecut=ecut):
                cprint("[%d]: ignoring ecut=%s because it's already in the DOJO_REPORT" % (dojo_trial, ecut), "magenta")
                continue

            # Build and register the work.
            pawecutdg = 2 * ecut if pseudo.ispaw else None
            work = phg_factory.work_for_pseudo(pseudo, kppa=1000, ecut=ecut, pawecutdg=pawecutdg,
                                               include_soc=options.soc)
            if work is not None:
                flow.register_work(work)
            else:
                cprint('Cannot create phgamma work for ecut %s, factory returned None' % str(ecut), "magenta")

    if len(flow) > 0:
        return flow.allocate()
    else:
        # Empty flow since all trials have been already performed.
        return None

@prof_main
def main():
    def str_examples():
        return """\
Usage Example:
    ppdojo_run.py Si.psp8  => Build pseudo_dojo flow for Si.fhi
"""

    def show_examples_and_exit(error_code=1):
        """Display the usage of the script."""
        print(str_examples())
        sys.exit(error_code)

    parser = argparse.ArgumentParser(epilog=str_examples())

    parser.add_argument('-m', '--manager', type=str, default=None,  help="Manager file")
    parser.add_argument('-d', '--dry-run', default=False, action="store_true",
                        help="Dry run, build the flow without submitting it")
    parser.add_argument('--paral-kgb', type=int, default=0,  help="Paral_kgb input variable.")
    parser.add_argument('-p', '--plot', default=False, action="store_true", help="Plot convergence when the flow is done")
    parser.add_argument('-n', '--new-ecut', type=int, default=None, action="store",
                        help="Extend the ecut grid with the new-ecut point")
    parser.add_argument('--soc', default=False, action="store_true", help=(
                        "Perform non-collinear run (nspinor==2, kptopt=3). Pseudo must have spin-orbit characteristic"))

    def parse_trials(s):
        if s == "all": return ["df", "gbrv", "phonon", "phowa"]
        return s.split(",")

    parser.add_argument('--trials', default="all",  type=parse_trials,
                        help=("List of tests e.g --trials=df,gbrv,phonon,phwoa\n"
                        "  df:     test delta factor against all electron reference.\n"
                        "  gbrv:   test fcc and bcc lattice parameters against AE reference.\n"
                        "  phonon: test convergence of phonon modes at gamma.\n"
                        "  phwoa:  test violation of the acoustic sum rule (without enforcing it) at the min and max ecut\n"))

    parser.add_argument('--loglevel', default="ERROR", type=str,
                        help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")
    parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
                         help='Verbose, can be supplied multiple times to increase verbosity')

    parser.add_argument('path', help='pseudopotential file.')

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

    options.manager = abilab.TaskManager.as_manager(options.manager)

    pseudo = dojopseudo_from_file(options.path)
    if pseudo is None:
        cprint("Error while parsing: %s" % options.path)
        cprint("Check your file. Returning 1")
        return 1

    flow = build_flow(pseudo, options)
    if flow is None:
        cprint("DOJO_REPORT is already computed for pseudo %s." % options.path, "magenta")
        return 1

    if options.dry_run:
        flow.build_and_pickle_dump()
        return 0
    else:
        # Run the flow with the scheduler.
        return flow.make_scheduler().start()


if __name__ == "__main__":
    sys.exit(main())
