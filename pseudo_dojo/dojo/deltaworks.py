from __future__ import division, print_function

import sys
import os
import collections
import numpy as np

from pymatgen.io.abinitio.pseudos import Pseudo
from pymatgen.io.abinitio.workflows import DeltaFactorWorkflow
from pseudo_dojo.refdata.deltafactor import df_database


class DeltaFactoryError(Exception):
    """Base Error class."""


class CIFNotFoundError(DeltaFactoryError):
    """CIF file not found in CIFs directory"""


class DeltaFactory(object):
    """
    Factory class producing work objects for the computation of the delta factor.
    """
    Error = DeltaFactoryError

    def __init__(self):
        self._dfdb = df_database()

    def get_cif_path(self, symbol):
        """Returns the path to the CIF file associated to the given symbol."""
        try:
            return self._dfdb.get_cif_path(symbol)
        except KeyError:
            raise CIFNotFoundError("%s: cannot find CIF file for symbol" % symbol)

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

        try:
            pseudo = Pseudo.aspseudo(pseudo)

            if pseudo.ispaw and pawecutdg is None:
                raise ValueError("pawecutdg must be specified for PAW calculations.")

            symbol = pseudo.symbol
        except AttributeError:
            print('error in parsing')
            # symbol = 'Ge'


        try:
            cif_path = self.get_cif_path(symbol)
        except Exception as exc:
            raise CIFNotFoundError(str(exc))
         # Include spin polarization for O, Cr and Mn (antiferromagnetic)
        # and Fe, Co, and Ni (ferromagnetic).
        spin_mode = "unpolarized"
        if symbol in ["Fe", "Co", "Ni"]: spin_mode = "polarized"
        if symbol in ["O", "Cr", "Mn"]: spin_mode = "afm"

        work = DeltaFactorWorkflow(cif_path, pseudo, kppa,
                         spin_mode=spin_mode, toldfe=toldfe, smearing=smearing, 
                         accuracy=accuracy, ecut=ecut, pawecutdg=pawecutdg, ecutsm=0.05, 
                         workdir=workdir, manager=manager, **kwargs
                        )
        return work

