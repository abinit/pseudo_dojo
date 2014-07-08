from __future__ import division, print_function

import sys
import os
#import collections
import numpy as np

from pymatgen.io.abinitio.pseudos import Pseudo
from pymatgen.io.abinitio.tasks import RelaxTask 
from pymatgen.io.abinitio.workflows import Workflow #, DeltaFactorWorkflow
from pseudo_dojo.refdata.gbrv import gbrv_database

import logging
logger = logging.getLogger(__name__)


class GBRVRelaxFactory(object):
    """
    Factory class producing `Workflow` objects for GBRV calculations.
    """
    def __init__(self):
        self._db = gbrv_database()

    def gbrv_relax_work(self, pseudo, struct_type, accuracy="normal", ngkpt=(8,8,8), 
        ecut=None, pawecutdg=None, toldfe=1.e-8, spin_mode="unpolarized", smearing="fermi_dirac:0.0005", 
        workdir=None, manager=None, **kwargs):
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

        # Get the entries in the database
        symbol = pseudo.symbol
        entry = self._db.get_entry(self, symbol, struct_type)

        # Build structures and handle missing values.
        structure = entry.build_structure(ref="ae")
        if structure is None:
            logger.warning("No AE structure for %s\n Will use gbrv_uspp data." % symbol)
            structure = fcc_entry.build_structure(ref="gbrv_uspp")

        assert structure is not None
        structure = AbiStructure.asabistructure(structure)

        return work



class GbrvRelaxWorkflow(Workflow):

    def __init__(self, structure, pseudo, ngkpt,
                 spin_mode="polarized", toldfe=1.e-8, smearing="fermi_dirac:0.1 eV",
                 accuracy="normal", ecut=None, pawecutdg=None, ecutsm=0.05, chksymbreak=0,
                 workdir=None, manager=None, **kwargs):
                 # FIXME Hack in chksymbreak
        """
        Build a `Workflow` for the computation of the deltafactor.

        Args:   
            structure:
                Structure object 
            pseudo:
                String with the name of the pseudopotential file or `Pseudo` object.
            nkkpt:
                MP divnsions.
            spin_mode:
                Spin polarization mode.
            toldfe:
                Tolerance on the energy (Ha)
            smearing:
                Smearing technique.
            workdir:
                String specifing the working directory.
            manager:
                `TaskManager` responsible for the submission of the tasks.
        """
        super(GbrvRelaxWorkflow, self).__init__(workdir=workdir, manager=manager)

        self.smearing = Smearing.assmearing(smearing)
        self.spin_mode = spin_mode
        self.accuracy = accuracy

        self.pseudo = Pseudo.aspseudo(pseudo)

        #extra_abivars = dict(
        #    pawecutdg=pawecutdg,
        #    ecutsm=ecutsm,
        #    toldfe=toldfe,
        #    prtwf=0,
        #    paral_kgb=0,
        #)
        #                               
        #extra_abivars.update(**kwargs)
                                                                                              
        #scf_input = ScfStrategy(new_structure, pseudo, ksampling,
        #                        accuracy=accuracy, spin_mode=spin_mode,
        #                        smearing=smearing, **extra_abivars)

        #                                                                                      
        #self.relax_task = RelaxTask(relax_input
        #work.register(scf_input, task_class=RelaxTask)

        #extra_abivars = dict(
        #    pawecutdg=pawecutdg,
        #    ecutsm=ecutsm,
        #    toldfe=toldfe,
        #    prtwf=0,
        #    paral_kgb=0,
        #)
        #                               
        #extra_abivars.update(**kwargs)
                                                                                                                          
        #self.ksampling = KSampling.monkhorst(ngkpt, chksymbreak=chksymbreak)
        #self.ksamping = Ksampling.gamma_centered(cls, kpts=nkkpt, use_symmetries=True, use_time_reversal=True):

        #relax_algo = 
        #relax_input = RelaxStrategy(structure, pseudos, ksampling, relax_algo, accuracy="normal", spin_mode="polarized",
        #                    smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None, **extra_abivars):
                                                                                                                          
        #scf_input = ScfStrategy(new_structure, pseudo, self.ksampling,
        #                        accuracy=accuracy, spin_mode=spin_mode,
        #                        smearing=smearing, **extra_abivars)


    def on_all_ok(self):
        # Get the relaxed structure.
        structure = self.relax_task.final_structure

        # Build tasks for the computation of lattice parameters and EOS.

        # From 94% to 106% of the equilibrium volume determined previously.
        self.volumes = structure.volume * np.arange(94, 108, 2) / 100.
                                                                                             
        for vol in self.volumes:
            new_lattice = structure.lattice.scale(vol)
                                                                                             
            new_structure = AbiStructure(new_lattice, structure.species, structure.frac_coords)

            #extra_abivars = dict(
            #    pawecutdg=pawecutdg,
            #    ecutsm=ecutsm,
            #    toldfe=toldfe,
            #    prtwf=0,
            #    paral_kgb=0,
            #)
            #                                                                                 
            #extra_abivars.update(**kwargs)
                                                                                             
            scf_input = ScfStrategy(new_structure, self.pseudo, self.ksampling,
                                    accuracy=self.accuracy, spin_mode=self.spin_mode,
                                    smearing=self.smearing, **extra_abivars)
                                                                                             
            # Register new task
            self.register(scf_input, task_class=ScfTask)

        return super(GbrvRelaxWorkflow, self).on_all_ok()
