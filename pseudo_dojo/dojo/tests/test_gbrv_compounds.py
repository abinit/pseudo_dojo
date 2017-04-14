"""Unit tests for dojo works"""
from __future__ import division, print_function, unicode_literals, absolute_import

import os
import pseudo_dojo.data as pdj_data

from abipy import abilab
from pseudo_dojo.core import PseudoDojoTest
from pseudo_dojo.dojo.gbrv_compounds import *


class GbrvCompoundTest(PseudoDojoTest):

    def test_nc_lif_gbrv_factory(self):
        """Testing GBRV work for SiO (rocksalt structure) with NC pseudos."""
        flow = abilab.Flow.temporary_flow()
        gbrv_factory = GbrvCompoundsFactory(xc="PBE")

        formula, struct_type, accuracy = "SiO", "rocksalt", "normal"
        pseudos = pdj_data.pseudos("Si.psp8", "O.psp8")

        work = gbrv_factory.relax_and_eos_work(accuracy, pseudos, formula, struct_type,
                                               ecut=6, pawecutdg=None)
        flow.register_work(work)
        flow.build_and_pickle_dump(abivalidate=True)
        print("Working in ", flow.workdir)
        assert len(flow[0]) == 1

        flow.rmtree()
