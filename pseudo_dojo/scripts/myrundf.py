#!/usr/bin/env python
"""Compute the deltafactor for a given pseudopotential."""
from __future__ import division, print_function

__author__ = 'setten'

import os
import sys
import abipy.abilab as abilab

from pseudo_dojo.dojo.works import DeltaFactory
from pymatgen.io.abinitio.pseudos import Pseudo
from pymatgen.core.periodic_table import PeriodicTable


def build_flow(pseudo, manager, accuracies=None):
    pseudo = Pseudo.from_file(pseudo)
    # Instantiate the TaskManager.

    factory = DeltaFactory()
    #extra = {}
    workdir = pseudo.basename + "_DFLOW"
    #if os.path.exists(workdir):
    #       raise ValueError("%s exists" % workdir)

    flow = abilab.Flow(workdir=workdir, manager=manager, pickle_protocol=0)
    kppa = 6750  # Use this to have the official k-point sampling

    extra_abivars = {
            "mem_test": 0,
            "fband": 2,
            "nstep": 100,
            "paral_kgb": 0,
            #"paral_kgb": 1,
            #"nsym": 1,
            #"nsppol": 2,
            #"nspden": 2,
    }

    report = pseudo.read_dojo_report()
    hints = report["hints"]
    #print(report)

    if accuracies is None:
        accuracies = ["normal", "high",]

    for accuracy in accuracies:
        ecut = hints[accuracy]["ecut"]
        pawecutdg = ecut * 2
        work = factory.work_for_pseudo(pseudo, kppa=kppa, ecut=ecut, pawecutdg=pawecutdg, toldfe=1.e-8, **extra_abivars)
        work.set_dojo_accuracy(accuracy)

        # Register the workflow.
        flow.register_work(work, workdir='W' + str(ecut))

    return flow.allocate()


def fireflow(flow):
    from fireworks import FireTaskBase, FWAction, Firework, LaunchPad, ScriptTask
    from fireworks.utilities.fw_serializers import FWSerializable
    from fireworks.core.rocket_launcher import launch_rocket
    from abipy.fworks.tasks import FireTaskWithFlow

    # set up the LaunchPad and reset it
    launchpad = LaunchPad()
    launchpad.reset('', require_password=False)

    # Build the flow
    #flow = build_flow()
    flow.build_and_pickle_dump()

    # create the Firework consisting of a single task
    firework = Firework(FireTaskWithFlow(flow=flow))

    # store workflow
    launchpad.add_wf(firework)

    #launch it locally
    #launch_rocket(launchpad)

    return 0

#@abilab.flow_main
def main():
    # lglevel is bound to the string value obtained from the command line argument.
    # Convert to upper case to allow the user to specify --loglevel=DEBUG or --loglevel=debug
    import logging
    #numeric_level = getattr(logging, options.loglevel.upper(), None)
    numeric_level = getattr(logging, "ERROR", None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % options.loglevel)
    logging.basicConfig(level=numeric_level)

    path = sys.argv[1]
    if len(sys.argv) >= 3:
        manager = abilab.TaskManager.from_file(sys.argv[2])
    else:
        manager = abilab.TaskManager.from_user_config()
        #print(manager)

    if os.path.isfile(path):
        #accuracies = ["normal",]
        accuracies = ["normal", "high",]
        accuracies = ["low", "normal", "high",]
        #accuracies = ["low",]
        flow = build_flow(path, manager, accuracies=accuracies)
        flow.build_and_pickle_dump()
        #flow.rapidfire()
        #print("nlaunch: %d" % flow.rapidfire())
        #fireflow(flow)
        #flow.make_scheduler().start()
    else:
        table = PeriodicTable()
        all_symbols = [element.symbol for element in table.all_elements]
        dirs = [os.path.join(path, d) for d in os.listdir(path) if d in all_symbols]
        #print("dirs", dirs)
        pseudos = []
        for d in dirs:
            print(d)
            pseudos.extend(os.path.join(d, p) for p in os.listdir(d) if p.endswith(".psp8"))
        print(pseudos)

        nflows, nlaunch = 0, 0
        for pseudo in pseudos:
            flow = build_flow(pseudo, manager, accuracies=["normal", "high",])
            if os.path.exists(flow.workdir) or nflows >= 4: continue
            nflows += 1
            flow.build_and_pickle_dump()
            nlaunch += flow.rapidfire()

        print("nlaunch: %d" % nlaunch)
        print("nflows: %d" % nflows)


if __name__ == "__main__":
    main()
