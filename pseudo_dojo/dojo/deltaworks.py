from __future__ import division, print_function

import sys
import os
import collections
import numpy as np

from pymatgen.io.abinitio.workflow import DeltaTest

##########################################################################################

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
        from pseudo_dojo.refdata.deltafactor import df_database
        self._dfdb = df_database()

    def get_cif_path(self, symbol):
        """Returns the path to the CIF file associated to the given symbol."""
        try:
            return self._dfdb.get_cif_path(symbol)
        except KeyError:
            raise CIFNotFoundError("%s: cannot find CIF file for symbol" % symbol)

    def work_for_pseudo(self, workdir, runmode, pseudo, accuracy="normal", kppa=6750, 
        ecut=None, toldfe=1.e-8, smearing="fermi_dirac:0.0005"):
        """
        Returns a Work object from the given pseudopotential.

        .. note: 
            0.001 Rydberg is the value used with WIEN2K
        """
        try:
            cif_path = self.get_cif_path(pseudo.symbol)
        except Exception as exc:
            raise CIFNotFoundError(str(exc))

        # Include spin polarization for O, Cr and Mn (antiferromagnetic) 
        # and Fe, Co, and Ni (ferromagnetic). 
        spin_mode = "unpolarized"

        if pseudo.symbol in ["Fe", "Co", "Ni"]: spin_mode = "polarized"
        if pseudo.symbol in ["O", "Cr", "Mn"]: spin_mode = "afm"

        work = DeltaTest(workdir, runmode, cif_path, pseudo, kppa,
                         spin_mode=spin_mode, toldfe=toldfe, smearing=smearing, 
                         accuracy=accuracy, ecut=ecut, ecutsm=0.05,
                        )
        return work

##########################################################################################
