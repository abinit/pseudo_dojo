"""Integration tests for deltafactor calculations with Abinit and AbiPy."""
from __future__ import print_function, division, unicode_literals, absolute_import

import pytest
import abipy.data as abidata
import abipy.abilab as abilab

from abipy.core.testing import has_abinit
from pseudo_dojo.dojo.works import DeltaFactory, GbrvFactory


def itest_deltafactor(fwp, tvars):
    """
    Testing the flow used for the computation of the deltafactor.
   """
    # Path of the pseudopotential to test.
    #pseudo = abidata.pseudo("Si.GGA_PBE-JTH-paw.xml").as_tmpfile()
    pseudo =  abilab.Pseudo.from_file("./Pb-d-3_r.psp8")

    flow = abilab.Flow(workdir=fwp.workdir, manager=fwp.manager)

    # Build the workflow for the computation of the deltafactor.
    # The workflow will produce a pdf file with the equation of state
    # and a file deltafactor.txt with the final results in the
    # outdir directory DELTAFACTOR/work_0/outdir.
    kppa = 20  # this value is for testing purpose (6570 is the correct one)
    ecut = 2
    pawecutdg = ecut * 2 if pseudo.ispaw else None

    work = DeltaFactory().work_for_pseudo(pseudo, kppa=kppa, ecut=ecut, pawecutdg=pawecutdg, 
                                          include_so=True,
                                          paral_kgb=tvars.paral_kgb)

    # Register the workflow.
    flow.register_work(work)
    flow.allocate()
    flow.build_and_pickle_dump()

    for task in flow[0]:
        task.start_and_wait()

    flow.check_status()
    flow.show_status()
    assert flow.all_ok
    assert all(work.finalized for work in flow)
    results = flow[0].get_results()
    #20.453 ang^3 88.545 GPa 4.31 20.8658081501 336.680999051 GPa -35.681897152
    #delta Equation of State: deltafactor_polyfit
    #Minimum volume = 20.87 Ang^3
    #modulus = 2.10 eV/Ang^3 = 336.68 GPa, b1 = -35.68
    #Deltafactor = 15.681 meV
    #assert 0
