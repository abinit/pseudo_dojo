"""Integration tests for pseudodojo."""
import pytest
import pseudo_dojo.data as pdj_data
import abipy.abilab as abilab

from pseudo_dojo.dojo.works import DeltaFactory, GbrvFactory


def itest_gbrv_gga_pawxml_flow(fwp, tvars):
    """Testing the GBRV flow with GGA and PAW-XML (relaxation + EOS)"""
    pseudo = pdj_data.pseudo("Si.GGA_PBE-JTH-paw.xml").as_tmpfile(tmpdir=fwp.workdir)
    assert pseudo is not None
    assert pseudo.has_dojo_report
    assert not pseudo.dojo_report.exceptions

    factory = GbrvFactory(pseudo.xc)

    ecut = 4
    pawecutdg = 2 * ecut if pseudo.ispaw else None
    flow = abilab.Flow(workdir=fwp.workdir, manager=fwp.manager)

    struct_types = ["fcc",] # "bcc"]
    assert not pseudo.dojo_report.has_trial("gbrv_fcc", ecut=ecut)

    for struct_type in struct_types:
        work = factory.relax_and_eos_work(pseudo, struct_type, ecut, pawecutdg=pawecutdg, paral_kgb=tvars.paral_kgb)
        flow.register_work(work)

    flow.build_and_pickle_dump(abivalidate=True)

    fwp.scheduler.add_flow(flow)
    assert fwp.scheduler.start() == 0
    assert not fwp.scheduler.exceptions

    flow.check_status(show=True)
    assert all(work.finalized for work in flow)
    assert flow.all_ok

    print(pseudo.dojo_report)
    assert pseudo.dojo_report.has_trial("gbrv_fcc", ecut=ecut)
    assert not pseudo.dojo_report.exceptions


def itest_gbrv_gga_ncsoc_flow(fwp, tvars):
    """Testing the GBRV flow with GGA and ONCVPSP+SO (relaxation + EOS)"""
    #return
    pseudo = pdj_data.pseudo("Pb-d-3_r.psp8").as_tmpfile(tmpdir=fwp.workdir)
    assert pseudo is not None
    assert pseudo.supports_soc
    assert not pseudo.dojo_report.exceptions

    factory = GbrvFactory(pseudo.xc)
    ecut = 2
    #ecut = 6
    pawecutdg = 2 * ecut if pseudo.ispaw else None

    flow = abilab.Flow(workdir=fwp.workdir, manager=fwp.manager)

    struct_types = ["fcc",] # "bcc"]
    assert pseudo.dojo_report.has_trial("gbrv_fcc", ecut=4)
    assert not pseudo.dojo_report.has_trial("gbrv_fcc", ecut=ecut)

    for struct_type in struct_types:
        work = factory.relax_and_eos_work(pseudo, struct_type, ecut, pawecutdg=pawecutdg,
                                          include_soc=True, paral_kgb=tvars.paral_kgb)
        flow.register_work(work)

    flow.build_and_pickle_dump(abivalidate=True)

    fwp.scheduler.add_flow(flow)
    assert fwp.scheduler.start() == 0
    assert not fwp.scheduler.exceptions

    flow.check_status(show=True)
    assert all(work.finalized for work in flow)
    assert flow.all_ok

    print(pseudo.dojo_report)
    assert pseudo.dojo_report.has_trial("gbrv_fcc_soc", ecut=ecut)
    assert not pseudo.dojo_report.exceptions
