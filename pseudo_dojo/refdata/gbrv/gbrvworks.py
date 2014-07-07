from __future__ import division, print_function

import sys
import os
import collections
import numpy as np

from pymatgen.io.abinitio.pseudos import Pseudo
from pymatgen.io.abinitio.workflows import DeltaFactorWorkflow
from pseudo_dojo.refdata.deltafactor import df_database


class GBRVFactoryError(Exception):
    """Base Error class."""


class GBRVFactory(object):
    """
    Factory class producing `Workflow` objects for GBRV calculations.
    """
    Error = GBRVFactoryError

    def __init__(self):
        self._gbrvdb = gbrv_database()

    def structure_for_pseudo(self, pseudo):
        pass

    def work_for_pseudo(self, pseudo, accuracy="normal", kppa=6750, 
        ecut=None, pawecutdg=None, toldfe=1.e-8, smearing="fermi_dirac:0.0005", 
        workdir=None, manager=None, **kwargs):
        """
        Returns a `Workflow` object from the given pseudopotential.

        Args:
            kwargs:
                Extra variables passed to Abinit.

        .. note: 
            0.001 Rydberg is the value used with WIEN2K
        """
        pseudo = Pseudo.aspseudo(pseudo)

        if pseudo.ispaw and pawecutdg is None:
            raise ValueError("pawecutdg must be specified for PAW calculations.")

        #try:
        #    cif_path = self.get_cif_path(pseudo.symbol)
        #except Exception as exc:
        #    raise CIFNotFoundError(str(exc))

        structure = self.structure_for_pseudo(pseudo)

        # Include spin polarization for O, Cr and Mn (antiferromagnetic) 
        # and Fe, Co, and Ni (ferromagnetic). 
        spin_mode = "unpolarized"

        #if pseudo.symbol in ["Fe", "Co", "Ni"]: spin_mode = "polarized"
        #if pseudo.symbol in ["O", "Cr", "Mn"]: spin_mode = "afm"

        #work = DeltaFactorWorkflow(cif_path, pseudo, kppa,
        #                 spin_mode=spin_mode, toldfe=toldfe, smearing=smearing, 
        #                 accuracy=accuracy, ecut=ecut, pawecutdg=pawecutdg, ecutsm=0.05, 
        #                 workdir=workdir, manager=manager, **kwargs)

        return work

    def on_all_ok(self):
        return self.get_results()

    def get_results(self):
        #wf_results = super(DeltaFactorWorkflow, self).get_results()


        wf_results.update({
            "etotal"    : list(etotal),
            "volumes"   : list(self.volumes),
            "natom"     : num_sites,
            "dojo_level": 1,
        })

        return wf_results
