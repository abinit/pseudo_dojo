"""Integration tests for EbandsFactory."""
from __future__ import print_function, division, unicode_literals, absolute_import

import pytest
import pseudo_dojo.data as pdj_data
import abipy.abilab as abilab

from pseudo_dojo.dojo.works import EbandsFactory


def itest_ebands_gga_pawxml_flow(fwp, tvars):
    """Testing the ebands flow for PAW-XML"""
    pseudo = pdj_data.pseudo("Si.GGA_PBE-JTH-paw.xml").as_tmpfile()
    assert pseudo is not None
    #assert not pseudo.has_dojo_report
    spin_mode = "unpolarized"

    flow = abilab.Flow(workdir=fwp.workdir, manager=fwp.manager)
    ebands_factory = EbandsFactory(xc=pseudo.xc)

    ecut = 8
    pawecutdg = 2 * ecut if pseudo.ispaw else None
    kppa = 20  # this value is for testing purpose 
    work = ebands_factory.work_for_pseudo(pseudo, accuracy="high", kppa=kppa, ecut=ecut, pawecutdg=pawecutdg, 
                                          spin_mode=spin_mode, bands_factor=2, smearing="fermi_dirac:0.0005",
                                          mem_test=0)
    flow.register_work(work)
    flow.build_and_pickle_dump()

    fwp.scheduler.add_flow(flow)
    assert fwp.scheduler.start() == 0
    assert not fwp.scheduler.exceptions

    flow.check_status(show=True)
    assert all(work.finalized for work in flow)
    assert flow.all_ok

    #assert pseudo.has_dojo_report
    #assert pseudo.dojo_report.has_trial("deltafactor")
    #assert pseudo.dojo_report.has_trial("deltafactor", ecut=ecut)


def itest_ebands_gga_ncsoc_flow(fwp, tvars):
    """Testing the ebands flow for NC+SOC pseudos."""
    #return
    pseudo = pdj_data.pseudo("Pb-d-3_r.psp8").as_tmpfile()
    assert pseudo is not None
    assert pseudo.supports_soc
    print(pseudo)
    spin_mode = "spinor"

    flow = abilab.Flow(workdir=fwp.workdir, manager=fwp.manager)
    ebands_factory = EbandsFactory(xc=pseudo.xc)

    ecut = 8
    pawecutdg = 2 * ecut if pseudo.ispaw else None
    kppa = 20  # this value is for testing purpose 
    work = ebands_factory.work_for_pseudo(pseudo, accuracy="high", kppa=kppa, ecut=ecut, pawecutdg=pawecutdg, 
                                          spin_mode=spin_mode, bands_factor=2, smearing="fermi_dirac:0.0005",
                                          mem_test=0)
    flow.register_work(work)
    flow.build_and_pickle_dump()

    fwp.scheduler.add_flow(flow)
    assert fwp.scheduler.start() == 0
    assert not fwp.scheduler.exceptions

    flow.check_status(show=True)
    assert all(work.finalized for work in flow)
    assert flow.all_ok

    #assert pseudo.has_dojo_report
    #assert pseudo.dojo_report.has_trial("deltafactor")
    #assert pseudo.dojo_report.has_trial("deltafactor", ecut=ecut)
