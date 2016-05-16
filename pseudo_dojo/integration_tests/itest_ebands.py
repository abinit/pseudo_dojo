"""Integration tests for pseudodojo."""
from __future__ import print_function, division, unicode_literals, absolute_import

import pytest
import abipy.data as abidata
import abipy.abilab as abilab

from pseudo_dojo.dojo.works import EbandsFactory

def itest_ebands(fwp, tvars):
    """Testing the ebands flow:"""
    pseudo = abidata.pseudo("Si.GGA_PBE-JTH-paw.xml").as_tmpfile()

    flow = abilab.Flow(workdir=fwp.workdir, manager=fwp.manager)
    ebands_factory = EbandsFactory(xc="PBE")

    ecut = 8
    pawecutdg = 2 * ecut if pseudo.ispaw else None
    kppa = 20  # this value is for testing purpose 
    work = ebands_factory.work_for_pseudo(pseudo, accuracy="high", kppa=kppa, ecut=ecut, pawecutdg=pawecutdg,
                                                  bands_factor=2, smearing="fermi_dirac:0.0005", qpt=[0,0,0], mem_test=0)

    flow.register_work(work)

    flow.allocate()
    flow.build_and_pickle_dump()

    fwp.scheduler.add_flow(flow)
    assert fwp.scheduler.start() == 0
    assert not fwp.scheduler.exceptions

    flow.check_status()
    flow.show_status()
    assert all(work.finalized for work in flow)
    assert flow.all_ok
