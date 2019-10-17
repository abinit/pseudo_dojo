"""Integration tests for GhostsFactory."""
import pytest
import pseudo_dojo.data as pdj_data
import abipy.abilab as abilab

from pseudo_dojo.dojo.works import GhostsFactory


def itest_ghosts_gga_pawxml_flow(fwp, tvars):
    """Testing the ghosts flow for PAW-XML"""
    pseudo = pdj_data.pseudo("Si.GGA_PBE-JTH-paw.xml").as_tmpfile(tmpdir=fwp.workdir)
    assert pseudo is not None
    print(pseudo)
    assert pseudo.has_dojo_report
    assert not pseudo.dojo_report.exceptions

    flow = abilab.Flow(workdir=fwp.workdir, manager=fwp.manager)
    factory = GhostsFactory(xc=pseudo.xc)

    ecut = 10
    pawecutdg = 2 * ecut if pseudo.ispaw else None
    # maxene 200 will make the task restart.
    work = factory.work_for_pseudo(pseudo, kppa=20, maxene=200, ecut=ecut, pawecutdg=pawecutdg,
                                   spin_mode="unpolarized")
    assert work.dojo_trial == "ghosts"
    flow.register_work(work)
    flow.build_and_pickle_dump(abivalidate=True)

    fwp.scheduler.add_flow(flow)
    assert fwp.scheduler.start() == 0
    assert not fwp.scheduler.exceptions

    flow.check_status(show=True)
    assert all(work.finalized for work in flow)
    assert flow.all_ok

    assert not pseudo.dojo_report.exceptions
    assert pseudo.dojo_report.has_trial("ghosts", ecut=ecut)
    key = "%.1f" % ecut
    data = pseudo.dojo_report["ghosts"][key]
    assert data["dojo_status"] == 0
    assert data["ecut"] == ecut
    ebands = abilab.ElectronBands.from_dict(data["ebands"])

    # This is to make sure that restart works as expected.
    #assert work[0].num_restarts > 0

    #fig = pseudo.dojo_report.plot_ebands(ecut=ecut, show=False)
    #assert fig is not None


def itest_ghosts_gga_ncsoc_flow(fwp, tvars):
    """Testing the ghosts flow for NC+SOC pseudos."""
    pseudo = pdj_data.pseudo("Pb-d-3_r.psp8").as_tmpfile(tmpdir=fwp.workdir)
    assert pseudo is not None
    print(pseudo)
    assert pseudo.has_dojo_report
    assert pseudo.supports_soc
    assert not pseudo.dojo_report.exceptions

    flow = abilab.Flow(workdir=fwp.workdir, manager=fwp.manager)
    factory = GhostsFactory(xc=pseudo.xc)

    ecut = 4
    pawecutdg = 2 * ecut if pseudo.ispaw else None
    kppa = 20  # this value is for testing purpose
    work = factory.work_for_pseudo(pseudo, kppa=kppa, maxene=2, ecut=ecut, pawecutdg=pawecutdg,
                                   spin_mode="spinor", include_soc=True)
    assert work.dojo_trial == "ghosts_soc"
    flow.register_work(work)
    flow.build_and_pickle_dump(abivalidate=True)

    fwp.scheduler.add_flow(flow)
    assert fwp.scheduler.start() == 0
    assert not fwp.scheduler.exceptions

    flow.check_status(show=True)
    assert all(work.finalized for work in flow)
    assert flow.all_ok

    # Reconstruct ElectronBands from JSON.
    assert not pseudo.dojo_report.exceptions
    print(list(pseudo.dojo_report.keys()))
    assert pseudo.dojo_report.has_trial("ghosts_soc", ecut=ecut)
    key = "%.1f" % ecut
    data = pseudo.dojo_report["ghosts_soc"][key]
    assert data["dojo_status"] == 0
    assert data["ecut"] == ecut
    ebands = abilab.ElectronBands.from_dict(data["ebands"])

    #fig = pseudo.dojo_report.plot_ebands(ecut=ecut, show=False)
    #assert fig is not None
