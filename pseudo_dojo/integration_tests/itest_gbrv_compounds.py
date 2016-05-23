"""Integration tests for pseudodojo."""
from __future__ import print_function, division, unicode_literals, absolute_import

import pytest
import pseudo_dojo.data as pdj_data
import abipy.abilab as abilab

from pseudo_dojo.dojo.gbrv_compounds import *


def itest_gbrvcompounds_gga_pawxml_flow(fwp, tvars):
    """Testing the GBRV flow with GGA and PAW-XML (relaxation + EOS)"""
    # Unit test
    flow = abilab.Flow(workdir=fwp.workdir, manager=fwp.manager)
    gbrv_factory = GbrvCompoundsFactory()

    extra_abivars = {
        "mem_test": 0,
        "fband": 2,
        "nstep": 100,
        "paral_kgb": 0,
    }

    formula, struct_type, accuracy = "LiF", "rocksalt", "normal"
    pseudos = pdj_data.pseudos("3li.pspnc", "9f.pspnc")

    work = gbrv_factory.relax_and_eos_work(accuracy, pseudos, formula, struct_type, 
                                           ecut=20, pawecutdg=None, **extra_abivars)
    flow.register_work(work)

    flow.build_and_pickle_dump()
    print("Working in ", flow.workdir)

    isok, results = flow.abivalidate_inputs()
    if not isok:
        print(results)
        assert isok
    assert len(flow[0]) == 1

    fwp.scheduler.add_flow(flow)
    assert fwp.scheduler.start() == 0
    assert not fwp.scheduler.exceptions
 
    flow.check_status()
    flow.show_status()
    assert all(work.finalized for work in flow)
    assert flow.all_ok
    assert len(flow[0]) == 1 + 9
