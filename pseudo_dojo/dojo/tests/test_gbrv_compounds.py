"""Unit tests for dojo works"""
from __future__ import division, print_function, unicode_literals

import os

from abipy import abilab
from abipy import data as abidata
from pseudo_dojo.core import PseudoDojoTest
from pseudo_dojo.dojo.gbrv_compounds import *


class GbrvCompoundTest(PseudoDojoTest):

    def test_nc_lif_gbrv_factory(self):
        """Testing GBRV work for NC LiF (rocksalt structure)."""
        flow = abilab.Flow.temporary_flow()
        gbrv_factory = GbrvCompoundsFactory()

        extra_abivars = {
            "mem_test": 0,
            "fband": 2,
            "nstep": 100,
            "paral_kgb": 0,
        }

        formula, struct_type, accuracy = "LiF", "rocksalt", "normal"
        pseudos = abidata.pseudos("3li.pspnc", "9f.pspnc")

        work = gbrv_factory.relax_and_eos_work(accuracy, pseudos, formula, struct_type, 
                                               ecut=6, pawecutdg=None, **extra_abivars)
        flow.register_work(work)

        flow.build_and_pickle_dump()
        print("Working in ", flow.workdir)

        isok, results = flow.abivalidate_inputs()
        if not isok:
            print(results)
            assert isok
        assert len(flow[0]) == 1

        flow.rmtree()


if __name__ == "__main__":
    import unittest
    unittest.main()
