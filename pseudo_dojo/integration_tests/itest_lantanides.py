"""Integration tests for GhostsFactory."""
import pytest
import pseudo_dojo.data as pdj_data
import abipy.abilab as abilab

from pseudo_dojo.dojo.works import RocksaltRelaxationFactory


def itest_lantanides_gga_flow(fwp, tvars):
    """Testing the ghosts flow for NC+SOC pseudos."""
    pseudo = pdj_data.pseudo("Lu-sp.psp8").as_tmpfile(tmpdir=fwp.workdir)
    assert pseudo is not None
    print(pseudo)
    assert pseudo.has_dojo_report
    assert not pseudo.supports_soc
    assert not pseudo.dojo_report.exceptions

    flow = abilab.Flow(workdir=fwp.workdir, manager=fwp.manager)
    factory = RocksaltRelaxationFactory(xc=pseudo.xc)

    ecut_list = [33]
    ngkpt = [1, 1, 1]
    #pawecutdg = 2 * ecut if pseudo.ispaw else None
    work = factory.work_for_pseudo(pseudo, ecut_list, pawecutdg=None, ngkpt=ngkpt, include_soc=False)
    assert work.dojo_trial == "raren_relax"
    flow.register_work(work)

    #flow.build_and_pickle_dump(abivalidate=True)
    assert flow.make_scheduler().start() == 0
    flow.check_status(show=True, verbose=1)

    assert all(work.finalized for work in flow)
    assert flow.all_ok

    # Reconstruct ElectronBands from JSON.
    #assert not pseudo.dojo_report.exceptions
    #print(list(pseudo.dojo_report.keys()))
    #assert pseudo.dojo_report.has_trial("ghosts_soc", ecut=ecut)
    #key = "%.1f" % ecut
    #data = pseudo.dojo_report["ghosts_soc"][key]
    #assert data["dojo_status"] == 0
    #assert data["ecut"] == ecut
    #ebands = abilab.ElectronBands.from_dict(data["ebands"])
