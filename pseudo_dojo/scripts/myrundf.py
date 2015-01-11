#!/usr/bin/env python
"""Compute the deltafactor for a given pseudopotential."""
from __future__ import division, print_function, unicode_literals

import sys
import os
import numpy as np
import abipy.abilab as abilab

from pseudo_dojo.dojo.works import DeltaFactory
from pymatgen.io.abinitio.pseudos import Pseudo
from pymatgen.core.periodic_table import PeriodicTable


def build_flow(pseudo, manager):
    pseudo = Pseudo.as_pseudo(pseudo)
    # Instantiate the TaskManager.

    factory = DeltaFactory()
    workdir = pseudo.basename + "_DFLOW"
    #if os.path.exists(workdir):
    #   raise ValueError("%s exists" % workdir)

    flow = abilab.Flow(workdir=workdir, manager=manager)
    # Use this to have the official k-point sampling
    kppa = 6750  

    extra_abivars = {
            "mem_test": 0,
            "fband": 2,
            "nstep": 100,
            #"paral_kgb": 0,
            "paral_kgb": 1,
            #"nsym": 1,
            #"nsppol": 2,
            #"nspden": 2,
    }

    report = pseudo.read_dojo_report()
    #print(report)
    #hints = report["hints"]
    ppgen_ecut = int(report["ppgen_hints"]["high"]["ecut"])

    #dense_right = np.linspace(ppgen_ecut, ppgen_ecut + 10, num=6)
    #dense_left = np.linspace(ppgen_ecut-8, ppgen_ecut, num=4, endpoint=False)
    #coarse_high = np.linspace(ppgen_ecut + 15, ppgen_ecut + 40, num=4)

    dense_right = np.arange(ppgen_ecut, ppgen_ecut + 6*2, step=2)
    dense_left = np.arange(max(ppgen_ecut-6, 2), ppgen_ecut, step=2)
    coarse_high = np.arange(ppgen_ecut + 15, ppgen_ecut + 35, step=5)

    ecut_list = list(dense_left) + list(dense_right) + list(coarse_high)

    for ecut in ecut_list:
        pawecutdg = 2 * ecut 
        # Build and register the workflow.
        work = factory.work_for_pseudo(pseudo, kppa=kppa, ecut=ecut, pawecutdg=pawecutdg, toldfe=1.e-8, **extra_abivars)
        flow.register_work(work, workdir='W' + str(ecut))

    return flow.allocate()


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
        flow = build_flow(path, manager)
        flow.build_and_pickle_dump()
        #flow.rapidfire()
        #print("nlaunch: %d" % flow.rapidfire())
        flow.make_scheduler().start()
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
            flow = build_flow(pseudo, manager)
            if os.path.exists(flow.workdir) or nflows >= 4: continue
            nflows += 1
            flow.build_and_pickle_dump()
            nlaunch += flow.rapidfire()

        print("nlaunch: %d" % nlaunch)
        print("nflows: %d" % nflows)

    return 0


if __name__ == "__main__":
    sys.exit(main())
