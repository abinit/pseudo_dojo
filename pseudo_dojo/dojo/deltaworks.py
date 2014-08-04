from __future__ import division, print_function

import numpy as np

from pymatgen.core.structure import Structure
from pymatgen.io.abinitio.pseudos import Pseudo
from pymatgen.io.smartio import read_structure
from pymatgen.io.gwwrapper.helpers import refine_structure
from pymatgen.io.abinitio.workflows import Workflow
from pymatgen.io.abinitio.abiobjects import AbiStructure, Smearing, KSampling, Electrons, RelaxationMethod
from pymatgen.io.abinitio.strategies import HtcStrategy, ScfStrategy, RelaxStrategy #, num_valence_electrons
from pymatgen.io.abinitio.tasks import (Task, AbinitTask, Dependency, Node, ScfTask, NscfTask, BseTask, RelaxTask)
from pymatgen.io.abinitio.eos import EOS
from pseudo_dojo.refdata.deltafactor import df_database, df_compute


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

    def work_for_pseudo(self, pseudo, accuracy="normal", kppa=6750, ecut=None, pawecutdg=None,
                        toldfe=1.e-8, smearing="fermi_dirac:0.0005", workdir=None, manager=None, **kwargs):
        """
        Returns a `Workflow` object from the given pseudopotential.

        Args:
            kwargs:
                Extra variables passed to Abinit.

        .. note: 
            0.001 Rydberg is the value used with WIEN2K
        """
        pseudo = Pseudo.as_pseudo(pseudo)
        symbol = pseudo.symbol

        if pseudo.ispaw and pawecutdg is None:
            raise ValueError("pawecutdg must be specified for PAW calculations.")

        try:
            cif_path = self.get_cif_path(symbol)
        except Exception as exc:
            raise CIFNotFoundError(str(exc))

        # Include spin polarization for O, Cr and Mn (antiferromagnetic)
        # and Fe, Co, and Ni (ferromagnetic).
        spin_mode = "unpolarized"
        if symbol in ["Fe", "Co", "Ni"]:
            spin_mode = "polarized"

        if symbol in ["O", "Cr", "Mn"]:
            spin_mode = "afm"
            if symbol == 'O':
                kwargs['spinat'] = [(0, 0, 1), (0, 0, -1)]
            elif symbol == 'Cr':
                kwargs['spinat'] = [(0, 0, 1), (0, 0, -1)]
            elif symbol == 'Mn':
                kwargs['spinat'] = [(0, 0, 1), (0, 0, -1), (0, 0, -1), (0, 0, 1)]

        return DeltaFactorWorkflow(
            cif_path, pseudo, kppa,
            spin_mode=spin_mode, toldfe=toldfe, smearing=smearing,
            accuracy=accuracy, ecut=ecut, pawecutdg=pawecutdg, ecutsm=0.05,
            workdir=workdir, manager=manager, **kwargs)


class DeltaFactorWorkflow(Workflow):
    """Workflow for the calculation of the deltafactor."""
    def __init__(self, structure_or_cif, pseudo, kppa,
                 spin_mode="polarized", toldfe=1.e-8, smearing="fermi_dirac:0.1 eV",
                 accuracy="normal", ecut=None, pawecutdg=None, ecutsm=0.05, chksymbreak=0,
                 paral_kgb=0, workdir=None, manager=None, **kwargs):
                 # FIXME Hack in chksymbreak
        """
        Build a `Workflow` for the computation of the deltafactor.

        Args:   
            structure_or_cif:
                Structure object or string with the path of the CIF file.
            pseudo:
                String with the name of the pseudopotential file or `Pseudo` object.`
            kppa:
                Number of k-points per atom.
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
        super(DeltaFactorWorkflow, self).__init__(workdir=workdir, manager=manager)

        if isinstance(structure_or_cif, Structure):
            structure = refine_structure(structure_or_cif, symprec=1e-6)
        else:
            # Assume CIF file
            structure = refine_structure(read_structure(structure_or_cif), symprec=1e-6)

        # Set extra_abivars
        extra_abivars = dict(
            pawecutdg=pawecutdg,
            ecutsm=ecutsm,
            toldfe=toldfe,
            prtwf=0,
            paral_kgb=paral_kgb,
        )

        extra_abivars.update(**kwargs)

        if ecut is not None:
            extra_abivars.update({"ecut": ecut})

        self.pseudo = Pseudo.as_pseudo(pseudo)

        structure = AbiStructure.asabistructure(structure)
        self._input_structure = structure
        v0 = structure.volume

        # From 94% to 106% of the equilibrium volume.
        self.volumes = v0 * np.arange(94, 108, 2) / 100.

        for vol in self.volumes:
            new_lattice = structure.lattice.scale(vol)

            new_structure = Structure(new_lattice, structure.species, structure.frac_coords)
            new_structure = AbiStructure.asabistructure(new_structure)

            ksampling = KSampling.automatic_density(new_structure, kppa,
                                                    chksymbreak=chksymbreak)

            scf_input = ScfStrategy(new_structure, self.pseudo, ksampling,
                                    accuracy=accuracy, spin_mode=spin_mode,
                                    smearing=smearing, **extra_abivars)

            self.register(scf_input, task_class=ScfTask, manager=manager)

    def get_results(self):
        wf_results = super(DeltaFactorWorkflow, self).get_results()

        num_sites = self._input_structure.num_sites
        etotals = self.read_etotals(unit="eV")

        wf_results.update(dict(
            etotals=list(etotals),
            volumes=list(self.volumes),
            num_sites=num_sites))

        try:
            #eos_fit = EOS.Murnaghan().fit(self.volumes/num_sites, etotals/num_sites)
            #print("murn",eos_fit)
            #eos_fit.plot(show=False, savefig=self.path_in_workdir("murn_eos.pdf"))

            # Use same fit as the one employed for the deltafactor.
            eos_fit = EOS.DeltaFactor().fit(self.volumes/num_sites, etotals/num_sites)
            #eos_fit.plot(show=False, savefig=self.outdir.path_in("eos.pdf"))

            # FIXME: This object should be moved to pseudo_dojo.
            # Get reference results (Wien2K).

            wien2k = df_database().get_entry(self.pseudo.symbol)
                                                                                                 
            # Compute deltafactor estimator.
            dfact = df_compute(wien2k.v0, wien2k.b0_GPa, wien2k.b1,
                               eos_fit.v0, eos_fit.b0_GPa, eos_fit.b1, b0_GPa=True)

            print("delta", eos_fit)
            print("Deltafactor = %.3f meV" % dfact)

            wf_results.update({
                "dfact_meV": dfact,
                "v0": eos_fit.v0,
                "b0": eos_fit.b0,
                "b0_GPa": eos_fit.b0_GPa,
                "b1": eos_fit.b1,
            })

        except EOS.Error as exc:
            wf_results.push_exceptions(exc)

        # Write data for the computation of the delta factor
        with open(self.outdir.path_in("deltadata.txt"), "w") as fh:
            fh.write("# Deltafactor = %s meV\n" % dfact)
            fh.write("# Volume/natom [Ang^3] Etotal/natom [eV]\n")
            for (v, e) in zip(self.volumes, etotals):
                fh.write("%s %s\n" % (v/num_sites, e/num_sites))

        return wf_results

    def on_all_ok(self):
        """Callback executed when all tasks in self have reached S_OK."""
        return self.get_results()
