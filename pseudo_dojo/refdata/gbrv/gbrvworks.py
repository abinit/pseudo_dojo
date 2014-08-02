from __future__ import division, print_function

from pymatgen.io.abinitio.pseudos import Pseudo
from pymatgen.io.abinitio.workflows import GbrvEosWorkflow, GbrvRelaxAndEosWorkflow
from pseudo_dojo.refdata.gbrv import gbrv_database

import abipy.abilab as abilab

import logging
logger = logging.getLogger(__name__)


class GbrvFactory(object):
    """
    Factory class producing `Workflow` objects for GBRV calculations.
    """
    def __init__(self):
        self._db = gbrv_database()

    def make_ref_structure(self, symbol, struct_type, ref):
        # Get the entry in the database
        entry = self._db.get_entry(symbol, struct_type)
                                                                                         
        # Build the structure and handle a possibly missing value.
        structure = entry.build_structure(ref=ref)

        if structure is None:
            logger.warning("No AE structure for %s\n Will use gbrv_uspp data." % symbol)
            structure = fcc_entry.build_structure(ref="gbrv_uspp")
        
        if structure is None: 
            logger.critical("Cannot initialize structure for %s, returning None!" % symbol)

        return structure

    def eoswork_for_pseudo(self, pseudo, struct_type, ecut, pawecutdg=None, paral_kgb=0, ref="ae"):
        pseudo = Pseudo.aspseudo(pseudo)
        if pseudo.ispaw and pawecutdg is None:
            raise ValueError("pawecutdg must be specified for PAW calculations.")

        structure = self.make_ref_structure(pseudo.symbol, struct_type=struct_type, ref=ref)

        return GbrvEosWorkflow(structure, struct_type, pseudo, ecut, pawecutdg=pawecutdg, paral_kgb=paral_kgb)

    def relax_and_eos_work(self, pseudo, struct_type, ecut, pawecutdg=None, paral_kgb=0, ref="ae"):
        """
        Returns a `Workflow` object from the given pseudopotential.

        Args:
            kwargs:
                Extra variables passed to Abinit.

        .. note: 
            GBRV tests are done with the following parameteres:

                - No spin polarization for structural relaxation 
                  (only for magnetic moments for which spin-unpolarized structures are used)
                - All calculations are done on an 8x8x8 k-point density and with 0.002 Ry Fermi-Dirac smearing
        """
        pseudo = Pseudo.aspseudo(pseudo)

        if pseudo.ispaw and pawecutdg is None:
            raise ValueError("pawecutdg must be specified for PAW calculations.")

        structure = self.make_ref_structure(pseudo.symbol, struct_type=struct_type, ref=ref)
 
        return GbrvRelaxAndEosWorkflow(structure, struct_type, pseudo, ecut,
                                       pawecutdg=pawecutdg, paral_kgb=paral_kgb)



def gbrv_flow_for_pseudo(pseudo, ecut, pawecutdg):
    manager = abilab.TaskManager.from_user_config()
    flow = abilab.AbinitFlow(workdir="FLOW_GBRV", manager=manager, pickle_protocol=0)

    struct_types = ["fcc"] #, "bcc"]
    #struct_types = ["bcc"]

    factory = GbrvFactory()

    for struct_type in struct_types:
        #work = factory.eoswork_for_pseudo(pseudo, struct_type, ecut, pawecutdg=pawecutdg) # ref="gbrv_paw")
        work = factory.relax_and_eos_work(pseudo, struct_type, ecut, pawecutdg=pawecutdg) #, ref="gbrv_paw")
        flow.register_work(work)

    flow.allocate()
    return flow


if __name__ == "__main__":
    pseudo = "si_pbe_v1_abinit.paw"
    pseudo = "o_pbe_v1.2_abinit.paw"
    pseudo = "c_pbe_v1.2_abinit.paw"
    pseudo = "ca_pbe_v1_abinit.paw"
    pseudo = "w_pbe_v1.2_abinit.paw"
    pseudo = "zn_pbe_v1_abinit.paw"

    import os
    pseudo = os.path.join("./GBRV_all_pbe", pseudo)
    ecut = 20
    pawecutdg = ecut * 4

    flow = gbrv_flow_for_pseudo(pseudo, ecut, pawecutdg)
    flow.build_and_pickle_dump()
