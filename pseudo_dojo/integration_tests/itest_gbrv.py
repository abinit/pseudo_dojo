"""Integration tests for pseudodojo."""
from __future__ import print_function, division, unicode_literals, absolute_import

import pytest
import abipy.data as abidata
import abipy.abilab as abilab

from abipy.core.testing import has_abinit
from pseudo_dojo.dojo.works import DeltaFactory, GbrvFactory

def itest_gbrv_flow(fwp, tvars):
    """Testing the GBRV flow: relaxation + EOS computation."""
    factory = GbrvFactory()

    #pseudo = "si_pbe_v1_abinit.paw"
    pseudo = abidata.pseudo("Si.GGA_PBE-JTH-paw.xml").as_tmpfile()
    ecut = 4
    pawecutdg = 2 * ecut if pseudo.ispaw else None

    flow = abilab.Flow(workdir=fwp.workdir, manager=fwp.manager)

    struct_types = ["fcc"] #, "bcc"]

    for struct_type in struct_types:
        work = factory.relax_and_eos_work(pseudo, struct_type, ecut, pawecutdg=pawecutdg, paral_kgb=tvars.paral_kgb)
        flow.register_work(work)

    flow.allocate()
    flow.build_and_pickle_dump()

    fwp.scheduler.add_flow(flow)
    assert fwp.scheduler.start() == 0
    assert not fwp.scheduler.exceptions

    #work = flow[0]
    #t0 = work[0]
    #assert len(work) == 1

    #t0.start_and_wait()
    #flow.check_status()

    # At this point on_all_ok is called.
    #assert t0.status == t0.S_OK
    #assert len(flow) == 2
    #assert len(flow[1]) == 9

    #assert not flow.all_ok

    #for task in flow[1]:
    #    task.start_and_wait()

    flow.check_status()
    flow.show_status()
    assert all(work.finalized for work in flow)
    assert flow.all_ok
