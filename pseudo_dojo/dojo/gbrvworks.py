from __future__ import division, print_function

import numpy as np
import abipy.abilab as abilab

from pymatgen.io.abinitio.eos import EOS
from pymatgen.io.abinitio.pseudos import Pseudo
from pymatgen.core.structure import Structure
from pymatgen.io.abinitio.abiobjects import AbiStructure, Smearing, KSampling, Electrons, RelaxationMethod
from pymatgen.io.abinitio.strategies import ScfStrategy, RelaxStrategy
from pymatgen.io.abinitio.tasks import ScfTask, RelaxTask
#from pymatgen.io.abinitio.workflows import Workflow
from pseudo_dojo.refdata.gbrv import gbrv_database
from pseudo_dojo.dojo.dojo_workflow import DojoWorkflow

import logging
logger = logging.getLogger(__name__)


__all__ = [
    "GbrvFactory",
    "GbrvRelaxAndEosWorkflow",
]


class GbrvFactory(object):
    """Factory class producing `Workflow` objects for GBRV calculations."""
    def __init__(self):
        self._db = gbrv_database()

    def make_ref_structure(self, symbol, struct_type, ref):
        """
        Return the structure used in the GBRV tests given the chemical symbol, the structure type
        and the reference code.
        """
        # Get the entry in the database
        entry = self._db.get_entry(symbol, struct_type)
                                                                                         
        # Build the structure and handle a possibly missing value.
        structure = entry.build_structure(ref=ref)

        if structure is None:
            logger.warning("No AE structure for %s\n Will use gbrv_uspp data." % symbol)
            structure = entry.build_structure(ref="gbrv_uspp")
        
        if structure is None: 
            logger.critical("Cannot initialize structure for %s, returning None!" % symbol)

        return structure

    def relax_and_eos_work(self, pseudo, struct_type, ecut=None, pawecutdg=None, paral_kgb=0, ref="ae"):
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
        pseudo = Pseudo.as_pseudo(pseudo)

        if pseudo.ispaw and pawecutdg is None:
            raise ValueError("pawecutdg must be specified for PAW calculations.")

        structure = self.make_ref_structure(pseudo.symbol, struct_type=struct_type, ref=ref)
 
        return GbrvRelaxAndEosWorkflow(structure, struct_type, pseudo, ecut=ecut,
                                       pawecutdg=pawecutdg, paral_kgb=paral_kgb)


def gbrv_nband(pseudo):
    # nband/fband are usually too small for the GBRV calculations.
    # FIXME this is not optimal
    nband = pseudo.Z_val
    nband += 0.5 * nband
    nband = int(nband)
    nband = max(nband,  8)
    print("nband", nband)
    return nband


class GbrvRelaxAndEosWorkflow(DojoWorkflow):

    def __init__(self, structure, struct_type, pseudo, ecut=None, pawecutdg=None, ngkpt=(8,8,8),
                 spin_mode="unpolarized", toldfe=1.e-8, smearing="fermi_dirac:0.001 Ha",
                 accuracy="normal", paral_kgb=0, ecutsm=0.05, chksymbreak=0,
                 workdir=None, manager=None, **kwargs):
                 # FIXME Hack in chksymbreak
        """
        Build a `Workflow` for the computation of the relaxed lattice parameter.

        Args:   
            structure:
                Structure object 
            structure_type:
                fcc, bcc 
            pseudo:
                String with the name of the pseudopotential file or `Pseudo` object.
            ecut:
                Cutoff energy in Hartree
            ngkpt:
                MP divisions.
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
        super(GbrvRelaxAndEosWorkflow, self).__init__(workdir=workdir, manager=manager)
        self.struct_type = struct_type
        self.accuracy = accuracy

        # nband must be large enough to accomodate fractional occupancies.
        self._pseudo = Pseudo.as_pseudo(pseudo)
        self.nband = gbrv_nband(self.pseudo)

        # Set extra_abivars.
        self.extra_abivars = dict(
            ecut=ecut,
            pawecutdg=pawecutdg,
            toldfe=toldfe,
            prtwf=0,
            #ecutsm=ecutsm,
            nband=self.nband,
            paral_kgb=paral_kgb)
                                       
        self.extra_abivars.update(**kwargs)
        self.ecut = ecut
        self.smearing = smearing

        self.ksampling = KSampling.monkhorst(ngkpt, chksymbreak=chksymbreak)
        self.spin_mode = spin_mode
        relax_algo = RelaxationMethod.atoms_and_cell()

        self.relax_input = RelaxStrategy(structure, pseudo, self.ksampling, relax_algo, 
                                         accuracy=accuracy, spin_mode=spin_mode,
                                         smearing=smearing, **self.extra_abivars)

        # Register structure relaxation task.
        self.relax_task = self.register(self.relax_input, task_class=RelaxTask)

    @property
    def dojo_trial(self):
        return "gbrv_" + self.struct_type

    @property
    def pseudo(self):
        return self._pseudo

    def add_eos_tasks(self):
        """
        Read the optimized structure from the netcdf file and add to self a new
        a new list of ScfTask for the computation of the EOS with the GBRV parameters.
        """
        # Get the relaxed structure.
        relaxed_structure = self.relax_task.read_final_structure()

        # GBRV use nine points from -1% to 1% of the initial guess and fitting the results to a parabola.
        # Note that it's not clear to me if they change the volume or the lattice parameter!
        self.volumes = relaxed_structure.volume * np.arange(99, 101.25, 0.25) / 100.

        for vol in self.volumes:
            new_lattice = relaxed_structure.lattice.scale(vol)
            new_structure = Structure(new_lattice, relaxed_structure.species, relaxed_structure.frac_coords)
            new_structure = AbiStructure.asabistructure(new_structure)

            scf_input = ScfStrategy(new_structure, self.pseudo, self.ksampling,
                                    accuracy=self.accuracy, spin_mode=self.spin_mode,
                                    smearing=self.smearing, **self.extra_abivars)

            # Register new task
            self.register(scf_input, task_class=ScfTask)

        # Allocate new tasks and update the pickle database.
        self.flow.allocate()
        self.flow.build_and_pickle_dump()

    def compute_eos(self):
        results = self.Results()

        # Read etotals and fit E(V) with a parabola to find minimum
        #num_sites = self._input_structure.num_sites
        etotals = self.read_etotals(unit="eV")[1:]
        assert len(etotals) == len(self.volumes)

        results.update(dict(
            etotals=list(etotals),
            volumes=list(self.volumes),
            #num_sites=num_sites,
        ))

        try:
            eos_fit = EOS.Quadratic().fit(self.volumes, etotals)
            #eos_fit.plot(show=False, savefig=self.outdir.path_in("eos.pdf"))

        except EOS.Error as exc:
            results.push_exceptions(exc)

        # Function to compute cubic a0 from primitive v0 (depends on struct_type)
        vol2a = {"fcc": lambda vol: (4 * vol) ** (1/3.),
                 "bcc": lambda vol: (2 * vol) ** (1/3.),
                 }[self.struct_type]

        a0 = vol2a(eos_fit.v0)

        results.update(dict(
            v0=eos_fit.v0,
            b0=eos_fit.b0,
            b1=eos_fit.b1,
            a0=a0,
            struct_type=self.struct_type,
        ))

        db = gbrv_database()
        entry = db.get_entry(self.pseudo.symbol, stype=self.struct_type)
        abs_err = a0 - entry.ae
        rel_err = 100 * (a0 - entry.ae) / entry.ae

        pawabs_err = a0 - entry.gbrv_paw
        pawrel_err = 100 * (a0 - entry.gbrv_paw) / entry.gbrv_paw

        print("for GBRV struct_type: ", self.struct_type, "a0= ", a0, "Angstrom")
        print("AE - THIS: abs_err = %f, rel_err = %f %%" % (abs_err, rel_err))
        print("GBRV-PAW - THIS: abs_err = %f, rel_err = %f %%" % (pawabs_err, pawrel_err))

        d = {k: results[k] for k in ("a0", "etotals", "volumes")}
        if results.exceptions:
            d["_exceptions"] = str(results.exceptions)

        self.write_dojo_report(d)

        return results

    @property
    def addeos_done(self):
        return len(self) > 1

    def on_all_ok(self):
        """
        This method is called when self reaches S_OK.
        It reads the optimized structure from the netcdf file and build
        a new workflow for the computation of the EOS with the GBRV parameters.
        """
        if not self.addeos_done:
            print("Building EOS tasks")
            self.add_eos_tasks()
            self._finalized = False
        else:
            print("Computing EOS")
            self.compute_eos()

        return super(GbrvRelaxAndEosWorkflow, self).on_all_ok()


def gbrv_flow_for_pseudo(workdir, pseudo, struct_type, ecut, pawecutdg, manager=None):

    manager = abilab.TaskManager.from_user_config() if manager is None else manager
    flow = abilab.AbinitFlow(workdir=workdir, manager=manager, pickle_protocol=0)

    factory = GbrvFactory()
    work = factory.relax_and_eos_work(pseudo, struct_type, ecut, pawecutdg=pawecutdg, ref="ae")

    flow.register_work(work)
    return flow.allocate()
