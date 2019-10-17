"""Integration tests for pseudodojo."""
import pytest
import pseudo_dojo.data as pdj_data
import abipy.abilab as abilab

from pseudo_dojo.dojo.gbrv_compounds import *


def itest_gbrvcompounds_gga_pawxml_flow(fwp, tvars):
    """Testing the GBRV flow with GGA and PAW-XML (relaxation + EOS)"""
    # Unit test
    flow = abilab.Flow(workdir=fwp.workdir, manager=fwp.manager)
    gbrv_factory = GbrvCompoundsFactory(xc="PBE")

    formula, struct_type, accuracy = "SiO", "rocksalt", "normal"
    pseudos = pdj_data.pseudos("Si.psp8", "O.psp8")

    work = gbrv_factory.relax_and_eos_work(accuracy, pseudos, formula, struct_type,
                                           ecut=12, pawecutdg=None)
    flow.register_work(work)

    flow.build_and_pickle_dump(abivalidate=True)
    print("Working in ", flow.workdir)

    assert len(flow[0]) == 1

    fwp.scheduler.add_flow(flow)
    assert fwp.scheduler.start() == 0
    assert not fwp.scheduler.exceptions

    flow.check_status(show=True)
    assert all(work.finalized for work in flow)
    assert flow.all_ok
    assert len(flow[0]) == 1 + 9
