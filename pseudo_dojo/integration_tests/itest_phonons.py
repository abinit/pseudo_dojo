"""Integration tests for pseudodojo."""
import pytest
import abipy.abilab as abilab
import pseudo_dojo.data as pdj_data

from pseudo_dojo.dojo.works import GammaPhononFactory


def itest_nc_phonons_gamma(fwp, tvars):
    """Testing the calculation of phonons at Gamma with/without the Asr (NC pseudos)."""
    #pseudo = pdj_data.pseudo("Si.GGA_PBE-JTH-paw.xml").as_tmpfile(fwp.workdir)
    pseudo = pdj_data.pseudo("Si.psp8").as_tmpfile(fwp.workdir)
    assert pseudo is not None
    assert pseudo.has_dojo_report
    assert not pseudo.dojo_report.exceptions

    flow = abilab.Flow(workdir=fwp.workdir, manager=fwp.manager)

    ecut = 8
    pawecutdg = 2 * ecut if pseudo.ispaw else None

    # This one requires deltafactor!
    factory = GammaPhononFactory(xc=pseudo.xc)
    work = factory.work_for_pseudo(pseudo, kppa=20, ecut=ecut, pawecutdg=pawecutdg)

    flow.register_work(work)
    flow.build_and_pickle_dump(abivalidate=True)

    #return
    fwp.scheduler.add_flow(flow)
    assert fwp.scheduler.start() == 0
    assert not fwp.scheduler.exceptions

    flow.check_status(show=True)
    assert all(work.finalized for work in flow)
    assert flow.all_ok

    assert not pseudo.dojo_report.exceptions
    assert pseudo.dojo_report.has_trial("phgamma", ecut=ecut)
