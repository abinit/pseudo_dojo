#!/usr/bin/env python
"""Compute the deltafactor for a given pseudopotential."""
from __future__ import division, print_function, unicode_literals

__author__ = 'setten'


import os
import sys
import abipy.abilab as abilab
from pseudo_dojo.dojo.dojo_workflows import DeltaFactorWorkflow
from pseudo_dojo.dojo.dojo_workflows import DeltaFactory


def build_flow(options):
    # Path of the pseudopotential to test.
    #pseudo = data.pseudo("14si.pspnc")
    #pseudo = data.pseudo("Si.GGA_PBE-JTH-paw.xml")
    here = os.path.abspath(os.path.curdir)
    ps_name = options['name']

    # the ocvpsps output file
    pseudo = os.path.join(here, ps_name+".out")
    print(pseudo)

    with open(pseudo, 'r') as fi:
        lines = fi.readlines()

    fo = open(ps_name+'.psp8', 'w')
    data = False
    for line in lines:
        if data:
            fo.write(line)
        if 'Begin PSPCODE8' in line:
            data = True
    fo.close()

    pseudo = os.path.join(here, ps_name+'.psp8')

    if options['strip']:
        sys.exit()

    # Instantiate the TaskManager.
    manager = abilab.TaskManager.from_user_config()  # if not options.manager else options.manager

    # Build the workflow for the computation of the deltafactor.
    # The calculation is done with the parameters and the cif files
    # used in the original paper. We only have to specify
    # the cutoff energy ecut (Ha) for the pseudopotential.
    # The workflow will produce a pdf file with the equation of state
    # and a file deltafactor.txt with the final results in the
    # outdir directory DELTAFACTOR/Wnn/outdir.

    factory = DeltaFactory()

    #extra = {}

    if options['test']:
        workdir = 'df_run_test_'+ps_name
        flow = abilab.AbinitFlow(workdir=workdir, manager=manager, pickle_protocol=0)
        kppa = 1000
        ecut = 40
        pawecutdg = ecut * 2
        work = factory.work_for_pseudo(pseudo, accuracy="normal", kppa=kppa,
                                       ecut=ecut, pawecutdg=pawecutdg,
                                       toldfe=1.e-8, smearing="fermi_dirac:0.0005")
        flow.register_work(work, workdir='W'+str(ecut))
    else:
        workdir = 'df_run_full_'+ps_name
        flow = abilab.AbinitFlow(workdir=workdir, manager=manager, pickle_protocol=0)
        kppa = 6750  # Use this to have the official k-point sampling
        for ecut in [12, 16, 20, 24, 28, 32, 36, 40, 44, 48]:
            pawecutdg = ecut * 2
            work = factory.work_for_pseudo(pseudo, accuracy="high", kppa=kppa,
                                           ecut=ecut, pawecutdg=pawecutdg,
                                           toldfe=1.e-10, smearing="fermi_dirac:0.0005")
            #this needs to be done in the loop over ecut ...
            flow.register_work(work, workdir='W'+str(ecut))

    flow.allocate()

    return flow.build_and_pickle_dump()


#abilab.flow_main
def main(options):
    print(options)
    build_flow(options)


if __name__ == "__main__":
    my_options = {'test': False, 'strip': False}

    name = sys.argv[1]
    my_options['name'] = name
    print(name, sys.argv)
    if not os.path.isfile(name+'.out'):
        print('No such file: %s.\nThe the first argument should be the name of the pseudo.' % name)
        raise RuntimeError

    for arg in sys.argv:
        my_options.update({arg: True})

    print(my_options)
    main(options=my_options)
