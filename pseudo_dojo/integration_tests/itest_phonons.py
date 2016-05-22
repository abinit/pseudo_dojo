"""Integration tests for pseudodojo."""
from __future__ import print_function, division, unicode_literals, absolute_import

import pytest
import abipy.data as abidata
import abipy.abilab as abilab

from abipy.core.testing import has_abinit
from pseudo_dojo.dojo.works import DFPTPhononFactory

def itest_phonons_without_asr(fwp, tvars):
    """Testing the calculation of phonons without Asr."""
    #pseudo = abidata.pseudo("Si.GGA_PBE-JTH-paw.xml").as_tmpfile()
    pseudo = abidata.pseudo("Si.oncvpsp").as_tmpfile()

    flow = abilab.Flow(workdir=fwp.workdir, manager=fwp.manager)

    ecut = 8
    pawecutdg = 2 * ecut if pseudo.ispaw else None
    kppa = 20  # this value is for testing purpose 

    # This one requires deltafactor!
    phonon_factory = DFPTPhononFactory(xc=pseudo.xc)
    work = phonon_factory.work_for_pseudo(pseudo, accuracy="high", kppa=kppa, ecut=ecut, pawecutdg=pawecutdg,
                                          tolwfr=1.e-20, smearing="fermi_dirac:0.0005", qpt=[0,0,0], rfasr=0)
    assert work is None
    return
    print(work)
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
