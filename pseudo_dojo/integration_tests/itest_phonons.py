"""Integration tests for pseudodojo."""
from __future__ import print_function, division, unicode_literals, absolute_import

import pytest
import abipy.abilab as abilab
import pseudo_dojo.data as pdj_data

from pseudo_dojo.dojo.works import DFPTPhononFactory


def itest_phonons_without_asr(fwp, tvars):
    """Testing the calculation of phonons without Asr."""
    #pseudo = pdj_data.pseudo("Si.GGA_PBE-JTH-paw.xml").as_tmpfile()
    pseudo = pdj_data.pseudo("Si.psp8").as_tmpfile()
    assert pseudo is not None
    assert pseudo.has_dojo_report
    assert not pseudo.dojo_report.exceptions

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
    flow.build_and_pickle_dump()

    fwp.scheduler.add_flow(flow)
    assert fwp.scheduler.start() == 0
    assert not fwp.scheduler.exceptions

    flow.check_status(show=True)
    assert all(work.finalized for work in flow)
    assert flow.all_ok

    #assert pseudo.dojo_report.has_trial("deltafactor", ecut=ecut)
    assert not pseudo.dojo_report.exceptions
