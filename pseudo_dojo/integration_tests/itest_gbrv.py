"""Integration tests for pseudodojo."""
from __future__ import print_function, division, unicode_literals, absolute_import

import pytest
import abipy.data as abidata
import abipy.abilab as abilab

from pseudo_dojo.dojo.works import DeltaFactory, GbrvFactory

def itest_gbrv_gga_paw_flow(fwp, tvars):
    """Testing the GBRV flow with GGA and PAW (relaxation + EOS)"""
    #return
    factory = GbrvFactory()

    pseudo = abidata.pseudo("Si.GGA_PBE-JTH-paw.xml").as_tmpfile()
    assert pseudo is not None
    ecut = 4
    pawecutdg = 2 * ecut if pseudo.ispaw else None

    flow = abilab.Flow(workdir=fwp.workdir, manager=fwp.manager)

    struct_types = ["fcc",] # "bcc"]

    for struct_type in struct_types:
        work = factory.relax_and_eos_work(pseudo, struct_type, ecut, pawecutdg=pawecutdg, paral_kgb=tvars.paral_kgb)
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

def itest_gbrv_gga_ncsoc_flow(fwp, tvars):
     """Testing the GBRV flow with GGA and ONCVPSP+SO (relaxation + EOS)"""
    #return
    factory = GbrvFactory()
 
    pseudo = abilab.Pseudo.from_file("./Pb-d-3_r.psp8")
    assert pseudo is not None
    assert pseudo.supports_soc
    ecut = 4
    pawecutdg = 2 * ecut if pseudo.ispaw else None
 
    flow = abilab.Flow(workdir=fwp.workdir, manager=fwp.manager)
 
    struct_types = ["fcc",] # "bcc"]
 
    for struct_type in struct_types:
        work = factory.relax_and_eos_work(pseudo, struct_type, ecut, pawecutdg=pawecutdg, 
                                          include_soc=True, paral_kgb=tvars.paral_kgb)
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
