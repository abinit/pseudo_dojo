"""Unit tests for dojo works"""
from __future__ import division, print_function, unicode_literals, absolute_import

import os

from abipy import abilab
from pseudo_dojo.core import PseudoDojoTest
import pseudo_dojo.data as pdj_data
from pseudo_dojo.dojo.works import *


class DeltaFactorTest(PseudoDojoTest):

    def test_nc_silicon_df(self):
        """Testing df factory for NC silicon."""
        pseudo = pdj_data.pseudo("Si.psp8")
        df_factory = DeltaFactory(xc=pseudo.xc)

        extra_abivars = {
            "mem_test": 0,
            "fband": 2,
            "nstep": 100,
            "paral_kgb": 0,
        }

        flow = abilab.Flow.temporary_flow()
        work = df_factory.work_for_pseudo(pseudo, kppa=1, ecut=2, pawecutdg=None, **extra_abivars)
        flow.register_work(work)

        flow.build_and_pickle_dump(abivalidate=True)
        flow.check_status(show=True)
        flow.rmtree()


class GbrvTest(PseudoDojoTest):

    def test_nc_silicon_gbrv_factory(self):
        """Testing GBRV work for NC silicon."""
        pseudo = pdj_data.pseudo("Si.psp8")
        gbrv_factory = GbrvFactory(xc=pseudo.xc)

        flow = abilab.Flow.temporary_flow()
        for struct_type in ("fcc", "bcc"):
            work = gbrv_factory.relax_and_eos_work(pseudo, struct_type, ecut=3, pawecutdg=None)
            flow.register_work(work)

        flow.build_and_pickle_dump(abivalidate=True)
        flow.check_status(show=True)
        assert len(flow[0]) == 1
        flow.rmtree()


class GhostsTest(PseudoDojoTest):

    def test_nc_silicon_ghosts_factory(self):
        """Testing GhostsWork for NC silicon."""
        pseudo = pdj_data.pseudo("Si.psp8")
        factory = GhostsFactory(xc=pseudo.xc)

        flow = abilab.Flow.temporary_flow()
        work = factory.work_for_pseudo(pseudo, ecut=3, pawecutdg=None)
        flow.register_work(work)

        flow.build_and_pickle_dump(abivalidate=True)
        flow.check_status(show=True)
        flow.rmtree()


class GammaPhononTest(PseudoDojoTest):

    def test_nc_silicon_gammaphonon_factory(self):
        """Testing GammaPhononFactory for NC silicon."""
        pseudo = pdj_data.pseudo("Si.psp8")
        phgamma_factory = GammaPhononFactory(xc=pseudo.xc)

        flow = abilab.Flow.temporary_flow()
        work = phgamma_factory.work_for_pseudo(pseudo, ecut=3, pawecutdg=None)
        flow.register_work(work)

        flow.build_and_pickle_dump(abivalidate=True)
        flow.check_status(show=True)
        flow.rmtree()
