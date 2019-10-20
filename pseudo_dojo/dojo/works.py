# coding: utf-8
"""Base class for Dojo Workkflows."""
import abc
import os
import logging
import numpy as np

from monty.io import FileLock
from pymatgen.core.xcfunc import XcFunc
from abipy.core.structure import Structure
from abipy.abio.factories import ion_ioncell_relax_input
from abipy.flowtk.abiobjects import SpinMode, Smearing, KSampling, RelaxationMethod
from abipy.flowtk.tasks import RelaxTask
from abipy.flowtk.works import Work, RelaxWork, PhononWork
from abipy import abilab
from pseudo_dojo.core.dojoreport import DojoReport, dojo_dfact_results, dojo_gbrv_results
from pseudo_dojo.refdata.gbrv import gbrv_database
from pseudo_dojo.refdata.deltafactor import df_database
from pseudo_dojo.refdata.lantanides.database import raren_database

logger = logging.getLogger(__name__)


class DojoWork(Work):
    __metaclass__ = abc.ABCMeta

    @abc.abstractproperty
    def dojo_pseudo(self):
        """:class:`Pseudo` object"""

    @abc.abstractproperty
    def dojo_trial(self):
        """String identifying the DOJO trial. Used to write results in the DOJO_REPORT."""

    @property
    def djrepo_path(self):
        """Path to the djrepo file."""
        root, ext = os.path.splitext(self.dojo_pseudo.filepath)
        return root + ".djrepo"

    def add_entry_to_dojoreport(self, entry, overwrite_data=False, pop_trial=False):
        """
        Write/update the DOJO_REPORT section of the pseudopotential.
        Important paramenters such as the name of the dojo_trial and the energy cutoff
        are provided by the sub-class.
        Client code is responsible for preparing the dictionary with the data.

        Args:
            entry: Dictionary with results.
            overwrite_data: If False, the routine raises an exception if this entry is
                already filled.
            pop_trial: True if the trial should be removed before adding the new entry.
        """
        #root, ext = os.path.splitext(self.dojo_pseudo.filepath)
        #djrepo = root + ".djrepo"
        djrepo = self.djrepo_path
        self.history.info("Writing dojreport data to %s" % djrepo)

        # Update file content with Filelock.
        with FileLock(djrepo):
            # Read report from file.
            file_report = DojoReport.from_file(djrepo)

            # Create new entry if not already there
            dojo_trial = self.dojo_trial

            if pop_trial:
                file_report.pop(dojo_trial, None)

            if dojo_trial not in file_report:
                file_report[dojo_trial] = {}

            # Convert float to string with 1 decimal digit.
            dojo_ecut = "%.1f" % self.ecut

            # Check that we are not going to overwrite data.
            if dojo_ecut in file_report[dojo_trial]:
                if not overwrite_data:
                    raise RuntimeError("dojo_ecut %s already exists in %s. Cannot overwrite data" %
                            (dojo_ecut, dojo_trial))
                else:
                    file_report[dojo_trial].pop(dojo_ecut)

            # Update file_report by adding the new entry and write new file
            file_report[dojo_trial][dojo_ecut] = entry

            # Write new dojo report and update the pseudo attribute
            file_report.json_write()
            self._pseudo.dojo_report = file_report


class FactoryError(Exception):
    """Base Error class raised by Factory objects."""


class GhostsFactory(object):
    """
    Factory producing work objects for the calculation of ebands for testing for ghosts.
    """
    Error = FactoryError

    def __init__(self, xc):
        """xc is the exchange-correlation functional e.g. PBE, PW."""
        # Get a reference to the deltafactor database. Used to get a structure
        self._dfdb = df_database(xc)

    def get_cif_path(self, symbol):
        """Returns the path to the CIF file associated to the given symbol."""
        try:
            return self._dfdb.get_cif_path(symbol)
        except KeyError:
            raise self.Error("%s: cannot find CIF file for symbol" % symbol)

    def work_for_pseudo(self, pseudo, kppa=3000, maxene=250, ecut=None, pawecutdg=None,
                        spin_mode="unpolarized", include_soc=False,
                        tolwfr=1.e-12, smearing="fermi_dirac:0.1 eV", workdir=None, manager=None, **kwargs):
        """
        Returns a :class:`Work` object from the given pseudopotential.

        Args:
            pseudo: :class:`Pseudo` object.
            kppa: Number of k-points per reciprocal atom.
            ecut: Cutoff energy in Hartree
            pawecutdg: Cutoff energy of the fine grid (PAW only)
            spin_mode: Spin polarization option
            tolwfr: Tolerance on the residuals.
            smearing: Smearing technique.
            workdir: Working directory.
            manager: :class:`TaskManager` object.
            kwargs: Extra variables passed to Abinit.
        """
        if pseudo.ispaw and pawecutdg is None:
            raise ValueError("pawecutdg must be specified for PAW calculations.")

        if pseudo.xc != self._dfdb.xc:
            raise ValueError(
                "Pseudo xc differs from the XC used to instantiate the factory\n"
                "Pseudo: %s, Database: %s" % (pseudo.xc, self._dfdb.xc))
        try:
            cif_path = self.get_cif_path(pseudo.symbol)
        except Exception as exc:
            raise self.Error(str(exc))

        # DO NOT CHANGE THE STRUCTURE REPORTED IN THE CIF FILE.
        structure = Structure.from_file(cif_path, primitive=False)

        return GhostsWork(
            structure, pseudo, kppa, maxene,
            spin_mode=spin_mode, include_soc=include_soc, tolwfr=tolwfr, smearing=smearing,
            ecut=ecut, pawecutdg=pawecutdg, ecutsm=0.5,
            workdir=workdir, manager=manager, **kwargs)


class GhostsWork(DojoWork):
    """Work for the calculation of the deltafactor."""

    def __init__(self, structure, pseudo, kppa, maxene,
                 ecut=None, pawecutdg=None, ecutsm=0.5,
                 spin_mode="unpolarized", include_soc=False, tolwfr=1.e-15, smearing="fermi_dirac:0.1 eV",
                 chksymbreak=0, workdir=None, manager=None, **kwargs):
        """
        Build a :class:`Work` for the computation of a bandstructure to check for ghosts.

        Args:
            structure: :class:`Structure` object
            pseudo: String with the name of the pseudopotential file or :class:`Pseudo` object.
            kppa: Number of k-points per reciprocal atom.
            maxene: 250 eV
            spin_mode: Spin polarization mode.
            include_soc=True of SOC should be included.
            tolwfr: Stopping criterion.
            smearing: Smearing technique.
            workdir: String specifing the working directory.
            manager: :class:`TaskManager` responsible for the submission of the tasks.
        """
        super(GhostsWork, self).__init__(workdir=workdir, manager=manager)
        self._pseudo = pseudo
        self.include_soc = include_soc

        spin_mode = SpinMode.as_spinmode(spin_mode)
        smearing = Smearing.as_smearing(smearing)

        # Here we find an initial guess for the number of bands
        # The goal is to reach maxene eV above the fermi level.
        # Assume ~ b4ev fact eV per band
        # Add a buffer of nbdbuf states and enforce an even number of states
        nval = structure.num_valence_electrons(self.dojo_pseudo)
        self.maxene = maxene
        b4ev = 1.3
        nband = int(nval + int(b4ev * self.maxene))
        #nband = nval // 2 + 10
        if spin_mode.nsppol == 1: nband // 2
        nbdbuf = max(int(0.1 * nband), 4)
        nband += nbdbuf
        nband += nband % 2

        # Set extra_abivars
        self.ecut, self.pawecutdg = ecut, pawecutdg

        extra_abivars = dict(
            ecut=ecut,
            pawecutdg=pawecutdg,
            ecutsm=ecutsm,
            nband=int(nband),
            nbdbuf=int(nbdbuf),
            tolwfr=tolwfr,
            #prtwf=0,
            nstep=200,
            chkprim=0,
            mem_test=0
        )

        extra_abivars.update(**kwargs)

        # Disable time-reversal if nspinor == 2
        self.kppa = kppa
        ksampling = KSampling.automatic_density(structure, kppa, chksymbreak=chksymbreak,
                                                use_time_reversal=spin_mode.nspinor == 1)

        scf_input = abilab.AbinitInput(structure=structure, pseudos=self.dojo_pseudo)
        scf_input.add_abiobjects(ksampling, smearing, spin_mode)
        scf_input.set_vars(extra_abivars)
        self.register_scf_task(scf_input)

        self.dojo_status = 0

    @property
    def dojo_pseudo(self):
        return self._pseudo

    @property
    def dojo_trial(self):
        if not self.include_soc:
            return "ghosts"
        else:
            return "ghosts_soc"

    def on_ok(self, sender):
        """
        This callback is called when one task reaches status S_OK.

        Here we extract the band structure from the GSR file and we save it in the JSON file.
        If maxene is not reached, task is set to unconverged so that the launcher will restart it
        and we can enter on_ok again.
        """
        task = self[0]
        # If this is not my business!
        if sender != task: return super(GhostsWork, self).on_ok(sender)

        # Read ebands from the GSR file, get also the minimum number of planewaves
        with task.open_gsr() as gsr:
            ebands = gsr.ebands
            min_npw = np.amin(gsr.reader.read_value("number_of_coefficients"))
            gsr_maxene = np.amax(ebands.eigens[:,:,-1] - ebands.fermie)
            gsr_nband = gsr.nband

        # Increase nband if we haven't reached maxene and restart
        # dojo_status:
        #     0: if run completed
        #    >0: convergence is not reached. Used when we save previous results before restarting.
        #    -1: Cannot increase bands anymore because we are close to min_npw.
        nband_sentinel = int(0.85 * min_npw)
        nband_sentinel += nband_sentinel % 2

        if gsr_maxene >= self.maxene:
            self.dojo_status = 0
            task.history.info("Convergence reached: gsr_maxene %s >= self.maxene %s" % (gsr_maxene, self.maxene))
        else:
            task.history.info("Convergence not reached. Will test if it's possible to restart the task.")
            task.history.info("gsr_maxene %s < self.maxene %s" % (gsr_maxene, self.maxene))
            nband = nband_old = int(task.input["nband"])

            if nband >= nband_sentinel:
                self.dojo_status = -1
                task.history.info("Reached maximum number of bands. Setting dojo_status to -1 and exit")
            else:
                # Use previous run to compute better estimate of b4ev.
                nval = task.input.num_valence_electrons
                b4ev = (nband_old - nval // 2) / gsr_maxene
                task.history.info("New b4ev: %s" % b4ev)
                nband = int(nval + int(b4ev * self.maxene))
                # Previous version.
                #nband += (0.2 * nband_old)
                nbdbuf = max(int(0.1 * nband), 4)
                nband += nbdbuf
                nband += nband % 2
                if nband >= nband_sentinel: nband = nband_sentinel
                self.dojo_status += 1

                # Restart.
                task.set_vars(nband=int(nband), nbdbuf=int(nbdbuf))
                task.restart()
                self.finalized = False

        # Convert to JSON and add results to the dojo report.
        entry = dict(ecut=self.ecut, pawecutdg=self.pawecutdg, min_npw=int(min_npw),
                     maxene_wanted=self.maxene, maxene_gsr=float(gsr_maxene), nband=gsr_nband,
                     dojo_status=self.dojo_status, kppa=self.kppa,
                     ebands=ebands.as_dict())

        # Use pop_trial to avoid multiple keys with the same value!
        self.add_entry_to_dojoreport(entry, pop_trial=True)

        return super(GhostsWork, self).on_ok(sender)


class DeltaFactory(object):
    """Factory class producing work objects for the computation of the delta factor."""
    Error = FactoryError

    def __init__(self, xc):
        """xc is the exchange-correlation functional e.g. PBE, PW."""
        # Get a reference to the deltafactor database
        self._dfdb = df_database(xc)

    def get_cif_path(self, symbol):
        """Returns the path to the CIF file associated to the given symbol."""
        try:
            return self._dfdb.get_cif_path(symbol)
        except KeyError:
            raise self.Error("%s: cannot find CIF file for symbol" % symbol)

    def work_for_pseudo(self, pseudo, kppa=6750, ecut=None, pawecutdg=None,
                        toldfe=1.e-9, smearing="fermi_dirac:0.1 eV", include_soc=False,
                        workdir=None, manager=None, **kwargs):
        """
        Returns a :class:`Work` object from the given pseudopotential.

        Args:
            pseudo: :class:`Pseudo` object.
            kppa: kpoint per reciprocal atom
            ecut: Cutoff energy in Hartree
            pawecutdg: Cutoff energy of the fine grid (PAW only)
            toldfe: Tolerance on the energy (Ha)
            smearing: Smearing technique.
            include_soc: True if pseudo has SO contributions and calculation should be done with nspinor=2
            workdir: String specifing the working directory.
            manager: :class:`TaskManager` responsible for the submission of the tasks.
            kwargs: Extra variables passed to Abinit.

        .. note::

            0.001 Rydberg is the value used with WIEN2K
        """
        symbol = pseudo.symbol
        if pseudo.ispaw and pawecutdg is None:
            raise ValueError("pawecutdg must be specified for PAW calculations.")

        if pseudo.xc != self._dfdb.xc:
            raise ValueError(
                "Pseudo xc differs from the XC used to instantiate the factory\n"
                "Pseudo: %s, Database: %s" % (pseudo.xc, self._dfdb.xc))

        try:
            cif_path = self.get_cif_path(symbol)
        except Exception as exc:
            raise self.Error(str(exc))

        # WARNING: DO NOT CHANGE THE STRUCTURE REPORTED IN THE CIF FILE.
        structure = Structure.from_file(cif_path, primitive=False)

        # Include spin polarization and initial spinat for particular elements
        # TODO: Not sure spinat is ok if LDA.
        kwargs["spinat"], spin_mode = self._dfdb.spinat_spinmode_for_symbol(symbol)
        if include_soc: spin_mode = "spinor"
        # This is needed for LDA
        if symbol in ("O", "Mn"):
            print("Got Oxygen or Mn")
            spin_mode = "polarized"

        # Magnetic elements:
        # Start from previous SCF run to avoid getting trapped in local minima
        connect = symbol in ("Fe", "Co", "Ni", "Cr", "Mn", "O", "Zn", "Cu")

        return DeltaFactorWork(
            structure, pseudo, kppa, connect,
            spin_mode=spin_mode, include_soc=include_soc, toldfe=toldfe, smearing=smearing,
            ecut=ecut, pawecutdg=pawecutdg, ecutsm=0.5,
            workdir=workdir, manager=manager, **kwargs)


class DeltaFactorWork(DojoWork):
    """Work for the calculation of the deltafactor."""

    def __init__(self, structure, pseudo, kppa, connect,
                 ecut=None, pawecutdg=None, ecutsm=0.5,
                 spin_mode="polarized", include_soc=False, toldfe=1.e-9, smearing="fermi_dirac:0.1 eV",
                 chksymbreak=0, workdir=None, manager=None, **kwargs):
        """
        Build a :class:`Work` for the computation of the deltafactor.

        Args:
            structure: :class:`Structure` object
            pseudo: :class:`Pseudo` object.
            kppa: Number of k-points per reciprocal atom.
            connect: True if the SCF run should be initialized from the previous run.
            spin_mode: Spin polarization mode.
            toldfe: Tolerance on the energy (Ha)
            smearing: Smearing technique.
            workdir: String specifing the working directory.
            manager: :class:`TaskManager` responsible for the submission of the tasks.
        """
        super(DeltaFactorWork, self).__init__(workdir=workdir, manager=manager)
        self._pseudo = pseudo
        self.include_soc = include_soc

        spin_mode = SpinMode.as_spinmode(spin_mode)
        smearing = Smearing.as_smearing(smearing)

        # Compute the number of bands from the pseudo and the spin-polarization.
        # Add 6 bands to account for smearing.
        #nval = structure.num_valence_electrons(self.pseudo)
        #nband = int(nval / spin_mode.nsppol) + 6

        # Set extra_abivars
        self.ecut, self.pawecutdg = ecut, pawecutdg

        extra_abivars = dict(
            ecut=ecut,
            pawecutdg=pawecutdg,
            ecutsm=ecutsm,
            toldfe=toldfe,
            prtwf=-1 if not connect else 1,
            chkprim=0,
            nstep=200,
            fband=2.0,   # 0.5 is the default value but it's not large enough from some systems.
            #paral_kgb=paral_kgb,
            #nband=nband,
            #mem_test=0,
        )

        extra_abivars.update(**kwargs)
        self._input_structure = structure
        v0 = structure.volume

        # From 94% to 106% of the equilibrium volume.
        self.volumes = v0 * np.arange(94, 108, 2) / 100.

        for vol in self.volumes:
            # Build new structure
            new_lattice = structure.lattice.scale(vol)
            new_structure = Structure(new_lattice, structure.species, structure.frac_coords)

            ksampling = KSampling.automatic_density(new_structure, kppa, chksymbreak=chksymbreak,
                                                    use_time_reversal=spin_mode.nspinor == 1)

            scf_input = abilab.AbinitInput(structure=new_structure, pseudos=self.dojo_pseudo)
            scf_input.add_abiobjects(ksampling, smearing, spin_mode)
            scf_input.set_vars(extra_abivars)

            # Magnetic materials with nspinor = 2 requires connection
            # and a double SCF run (nsppol = 2 first then nspinor = 2).
            if connect and spin_mode.nspinor == 2:
                print("Using collinear then noncollinear scf task")
                self.register_collinear_then_noncollinear_scf_task(scf_input)
            else:
                self.register_scf_task(scf_input)

        if connect:
            logger.info("Connecting SCF tasks using previous WFK file")
            middle = len(self.volumes) // 2
            filetype = "WFK"
            for i, task in enumerate(self[:middle]):
                task.add_deps({self[i + 1]: filetype})

            for i, task in enumerate(self[middle+1:]):
                task.add_deps({self[middle + i]: filetype})

    @property
    def dojo_pseudo(self):
        return self._pseudo

    @property
    def dojo_trial(self):
        if not self.include_soc:
            return "deltafactor"
        else:
            return "deltafactor_soc"

    def on_all_ok(self):
        """
        This method is called once the `Work` is completed i.e. when all the tasks
        have reached status S_OK. Here we gather the results of the different tasks,
        the deltafactor is computed and the results are stored in the JSON file.
        """
        num_sites = self._input_structure.num_sites
        etotals = self.read_etotals(unit="eV")

        d, eos_fit = dojo_dfact_results(self.dojo_pseudo, num_sites, self.volumes, etotals)

        print("[%s]" % self.dojo_pseudo.symbol, "eos_fit:", eos_fit)
        print("Ecut %.1f, dfact = %.3f meV, dfactprime %.3f meV" % (self.ecut, d["dfact_meV"], d["dfactprime_meV"]))

        self.add_entry_to_dojoreport(d)

        return dict(returncode=0, message="Delta factor computed")


class GbrvFactory(object):
    """Factory class producing :class:`Work` objects for GBRV calculations."""

    def __init__(self, xc):
        """xc: exchange-correlation functional e.g. PBE or PW."""
        self.db = gbrv_database(xc)

    def make_ref_structure(self, symbol, struct_type, ref):
        """
        Return the structure used in the GBRV tests given the chemical symbol,
        the structure type and the reference code.
        """
        # Get the entry in the database
        entry = self.db.get_entry(symbol, struct_type)

        if entry is None:
            logger.critical("Cannot find entry for %s, returning None!" % symbol)
            return None

        # Build the structure and handle a possibly missing value.
        structure = entry.build_structure(ref=ref)

        if structure is None:
            logger.warning("No AE structure for %s\n Will use gbrv_uspp data." % symbol)
            structure = entry.build_structure(ref="gbrv_uspp")

        if structure is None:
            logger.critical("Cannot initialize structure for %s, returning None!" % symbol)

        return structure

    @property
    def xc(self):
        return self.db.xc

    def relax_and_eos_work(self, pseudo, struct_type, ecut=None, pawecutdg=None, include_soc=False,
                           ref="ae", **kwargs):
        """
        Returns a :class:`Work` object from the given pseudopotential.

        Args:
            include_soc: True if pseudo has SO contributions and calculation should be done with nspinor=2
            kwargs: Extra variables passed to Abinit.

        .. note::

            GBRV tests are done with the following parameteres:

                - No spin polarization for structural relaxation
                  (only for magnetic moments for which spin-unpolarized structures are used)
                - All calculations are done on an 8x8x8 k-point density and with 0.002 Ry Fermi-Dirac smearing
        """
        if pseudo.ispaw and pawecutdg is None:
            raise ValueError("pawecutdg must be specified for PAW calculations.")

        if pseudo.xc != self.xc:
            raise ValueError(
                "Pseudo xc differs from the XC used to instantiate the factory\n"
                "Pseudo: %s, Database: %s" % (pseudo.xc, self.xc))

        # Select spin_mode from include_soc.
        spin_mode = "unpolarized"
        if include_soc:
            spin_mode = "spinor"
            if not pseudo.supports_soc:
                raise ValueError("Pseudo %s does not support SOC calculation." % pseudo)

        structure = self.make_ref_structure(pseudo.symbol, struct_type=struct_type, ref=ref)

        return GbrvRelaxAndEosWork(structure, struct_type, pseudo,
                                   ecut=ecut, pawecutdg=pawecutdg, spin_mode=spin_mode, include_soc=include_soc,
                                   **kwargs)


class GbrvRelaxAndEosWork(DojoWork):

    def __init__(self, structure, struct_type, pseudo, ecut=None, pawecutdg=None, ngkpt=(8, 8, 8),
                 spin_mode="unpolarized", include_soc=False, toldfe=1.e-9, smearing="fermi_dirac:0.001 Ha",
                 ecutsm=0.05, chksymbreak=0,
                 workdir=None, manager=None, **kwargs):
        """
        Build a :class:`Work` for the computation of the relaxed lattice parameter.

        Args:
            structure: :class:`Structure` object
            structure_type: fcc, bcc
            pseudo: :class:`Pseudo` object.
            ecut: Cutoff energy in Hartree
            ngkpt: MP divisions.
            spin_mode: Spin polarization mode.
            toldfe: Tolerance on the energy (Ha)
            smearing: Smearing technique.
            workdir: String specifing the working directory.
            manager: :class:`TaskManager` responsible for the submission of the tasks.
        """
        super(GbrvRelaxAndEosWork, self).__init__(workdir=workdir, manager=manager)
        self.struct_type = struct_type

        # nband must be large enough to accomodate fractional occupancies.
        self._pseudo = pseudo
        self.include_soc = include_soc

        def gbrv_nband(pseudo):
            # nband/fband are usually too small for the GBRV calculations.
            # FIXME this is not optimal
            nband = pseudo.Z_val
            nband += 0.5 * nband
            nband = int(nband)
            nband = max(nband,  8)
            # Use even numer of bands. Needed when nspinor == 2
            if nband % 2 != 0: nband += 1
            return nband

        nband = gbrv_nband(self.dojo_pseudo)

        # Set extra_abivars.
        self.extra_abivars = dict(
            ecut=ecut,
            pawecutdg=pawecutdg,
            toldfe=toldfe,
            prtwf=-1,
            nband=nband,
            #ecutsm=0.5,
            #paral_kgb=paral_kgb
        )

        self.extra_abivars.update(**kwargs)
        self.ecut = ecut
        self.smearing = Smearing.as_smearing(smearing)

        # Kpoint sampling: shiftk depends on struct_type
        shiftk = {"fcc": [0, 0, 0], "bcc": [0.5, 0.5, 0.5]}.get(struct_type)
        self.spin_mode = SpinMode.as_spinmode(spin_mode)
        self.ksampling = KSampling.monkhorst(ngkpt, chksymbreak=chksymbreak, shiftk=shiftk,
                                            use_time_reversal=self.spin_mode.nspinor == 1)
        relax_algo = RelaxationMethod.atoms_and_cell()

        inp = abilab.AbinitInput(structure, pseudo)
        inp.add_abiobjects(self.ksampling, relax_algo, self.spin_mode, self.smearing)
        inp.set_vars(self.extra_abivars)

        # Register structure relaxation task.
        self.relax_task = self.register_relax_task(inp)

    @property
    def dojo_trial(self):
        if not self.include_soc:
            return "gbrv_" + self.struct_type
        else:
            return "gbrv_" + self.struct_type + "_soc"

    @property
    def dojo_pseudo(self):
        return self._pseudo

    def add_eos_tasks(self):
        """
        Read the optimized structure from the netcdf file and add to self a new
        a new list of ScfTask for the computation of the EOS with the GBRV parameters.
        """
        self.history.info("Building EOS tasks")

        # Get the relaxed structure.
        self.relaxed_structure = relaxed_structure = self.relax_task.get_final_structure()

        # GBRV use nine points from -1% to 1% of the initial guess and fitting the results to a parabola.
        # Note that it's not clear to me if they change the volume or the lattice parameter!
        self.volumes = relaxed_structure.volume * np.arange(99, 101.25, 0.25) / 100.

        for vol in self.volumes:
            new_lattice = relaxed_structure.lattice.scale(vol)
            new_structure = Structure(new_lattice, relaxed_structure.species, relaxed_structure.frac_coords)

            # Add ecutsm
            extra = self.extra_abivars.copy()
            extra["ecutsm"] = 0.5

            scf_input = abilab.AbinitInput(new_structure, self.dojo_pseudo)
            scf_input.add_abiobjects(self.ksampling, self.spin_mode, self.smearing)
            scf_input.set_vars(extra)

            # Register new task
            self.register_scf_task(scf_input)

        # Allocate new tasks and update the pickle database.
        self.flow.allocate()
        self.flow.build_and_pickle_dump()

    def compute_eos(self):
        """
        This method compute the equation of state following the approach described in the GBRV paper.
        Build dictionary with results and insert in the DOJO REPORT.
        Return: dictionary with results
        """
        self.history.info("Computing EOS")

        # Read etotals and fit E(V) with a parabola to find the minimum
        etotals = self.read_etotals(unit="eV")[1:]
        assert len(etotals) == len(self.volumes)

        num_sites = len(self.relaxed_structure)
        dojo_entry, eos_fit = dojo_gbrv_results(self.dojo_pseudo, self.struct_type, num_sites, self.volumes, etotals)
        self.add_entry_to_dojoreport(dojo_entry)

    @property
    def add_eos_done(self):
        """True if the EOS has been computed."""
        return len(self) > 1

    def on_all_ok(self):
        """
        This method is called when self reaches S_OK. It reads the optimized structure
        from the netcdf file and builds a new work for the computation of the EOS
        with the GBRV parameters.
        """
        if not self.add_eos_done:
            # Build SCF tasks for the EOS and tell the world we are not done!
            self.add_eos_tasks()
            self.finalized = False
        else:
            # Compute EOS, write data and enter in finalized mode.
            self.compute_eos()

        return super(GbrvRelaxAndEosWork, self).on_all_ok()


class GammaPhononFactory(object):
    """
    Factory class producing a specialized `Workflow` for Relaxation + Phonon calculations at the Gamma point.
    The initial structural parameters are taken from the deltafactor database.

    The work returned by `work_for_pseudo` peforms a full structural relaxation: ions + ions_cell
    Once the optimized structure is known, a new Workflow for GS + DFPT with the optimized parameters
    is created and added to the flow. This second work merges the DDB files, calls anaddb to compute
    the acoustic frequencies with asr in [0, 2] and save the results in the dojoreport.

    Client code usually builds severals works with different values of ecut in order
    to monitor the convergence of the frequencies and the fulfillment of the accoustic sum rule.
    """
    Error = FactoryError

    def __init__(self, xc):
        """xc is the exchange-correlation functional e.g. PBE, PW."""
        # Get reference to the deltafactor database
        self._dfdb = df_database(xc)

    def get_cif_path(self, symbol):
        """Returns the path to the CIF file associated to the given symbol."""
        try:
            return self._dfdb.get_cif_path(symbol)
        except KeyError:
            raise self.Error("%s: cannot find CIF file for symbol" % symbol)

    def work_for_pseudo(self, pseudo, kppa=1000, ecut=None, pawecutdg=None,
                        smearing="fermi_dirac:0.1 eV", include_soc=False,
                        workdir=None, manager=None):
        """
        Create and return a :class:`RelaxAndAddPhGammaWork` object.

        Args:
            pseudo: filepath or :class:`Pseudo` object.
            kppa: Number of k-points per reciprocal atom.
            ecut: Cutoff energy in Hartree
            pawecutdg: Cutoff energy of the fine grid (PAW only)
            smearing: Smearing technique.
            include_soc=True of SOC should be included.
            workdir: Working directory.
            manager: :class:`TaskManager` object.
        """
        symbol = pseudo.symbol
        if pseudo.ispaw and pawecutdg is None:
            raise ValueError("pawecutdg must be specified for PAW calculations.")

        if pseudo.xc != self._dfdb.xc:
            raise ValueError(
                "Pseudo xc differs from the XC used to instantiate the factory\n"
                "Pseudo: %s, Database: %s" % (pseudo.xc, self._dfdb.xc))

        qpt = np.zeros(3)
        self.include_soc = include_soc

        try:
            cif_path = self.get_cif_path(pseudo.symbol)
        except Exception as exc:
            raise self.Error(str(exc))

        # DO NOT CHANGE THE STRUCTURE REPORTED IN THE CIF FILE.
        structure = Structure.from_file(cif_path, primitive=False)

        # Get spinat and spin_mode from df database.
        spinat, spin_mode = self._dfdb.spinat_spinmode_for_symbol(symbol)

        # DFPT with AFM is not supported. Could try AFM in Relax and then polarized
        # in WFK + DFPT but this one is safer.
        if spin_mode == "afm": spin_mode = "polarized"
        if include_soc: spin_mode = "spinor"

        # Build inputs for structural relaxation.
        multi = ion_ioncell_relax_input(
                            structure, pseudo,
                            kppa=kppa, nband=None,
                            ecut=ecut, pawecutdg=pawecutdg, accuracy="normal", spin_mode=spin_mode,
                            smearing=smearing)

        # Set spinat from internal database.
        multi.set_vars(chkprim=0, mem_test=0, spinat=spinat)

        # Construct a *specialized" work for structural relaxation
        # This work will create a new workflow for phonon calculations
        # with the final relaxed structure (see on_all_ok below).
        work = RelaxAndAddPhGammaWork(ion_input=multi[0], ioncell_input=multi[1])

        # Monkey patch work
        work.dojo_kppa = kppa
        work.dojo_qpt = qpt
        work.ecut = ecut
        work.dojo_pawecutdg = pawecutdg
        work.dojo_include_soc = include_soc
        work._dojo_trial = "phgamma" if not include_soc else "phgamma_soc"
        work.dojo_pseudo = pseudo

        return work


class RelaxAndAddPhGammaWork(RelaxWork):
    """
    Work for structural relaxations. The first task relaxes the atomic position
    while keeping the unit cell parameters fixed. The second task uses the final
    structure to perform a structural relaxation in which both the atomic positions
    and the lattice parameters are optimized.
    """

    def on_all_ok(self):
        """
        Here I extend the implementation of super in order to create a new workflow
        for phonons with the optimized structural parameters.
        """
        results = super(RelaxAndAddPhGammaWork, self).on_all_ok()

        # Get the relaxed structure.
        relax_task = self[1]
        final_structure = relax_task.get_final_structure()

        # Use new structure in GS + DFPT runs and change some values.
        scf_input = relax_task.input.deepcopy()
        scf_input.set_structure(final_structure)

        # Remove input variables that can enter into conflict with DFPT.
        scf_input.pop_tolerances()
        scf_input.pop_par_vars()
        scf_input.pop_irdvars()
        scf_input.pop_vars(["ionmov", "optcell", "ntime", "dilatmx"])
        scf_input.set_vars(tolwfr=1e-20, nstep=80, nbdbuf=4)
        #nval = scf_input.num_valence_electrons

        # Build GS work and Phonon Work
        work = PhononDojoWork.from_scf_input(scf_input, self.dojo_qpt)
        for task in work[1:]:
            task.set_vars(prtwf=-1)

        # Monkey patch work
        work.dojo_kppa = self.dojo_kppa
        work.dojo_qpt = self.dojo_qpt
        work.ecut = self.ecut
        work.dojo_pawecutdg = self.dojo_pawecutdg
        work.dojo_include_soc = self.dojo_include_soc
        work._dojo_trial = "phgamma" if not self.dojo_include_soc else "phgamma_soc"
        work._pseudo = self.dojo_pseudo

        self.flow.register_work(work)
        # Allocate new tasks and update the pickle database.
        self.flow.allocate()
        self.flow.build_and_pickle_dump()

        return results


class PhononDojoWork(PhononWork, DojoWork):

    @property
    def dojo_trial(self):
        return self._dojo_trial

    @property
    def dojo_pseudo(self):
        return self._pseudo

    def on_all_ok(self):
        """
        merge the DDB files, invoke anaddb to compute the phonon frequencies with/without ASR
        Results are written to the dojoreport.
        """
        # Call super to merge the DDB files.
        results = super(PhononDojoWork, self).on_all_ok()

        # Read the final DDB produced in outdir.
        out_ddb = self.outdir.path_in("out_DDB")
        ddb = abilab.DdbFile(out_ddb)

        # Call anaddb with/without ASR.
        asr2_phbands = ddb.anaget_phmodes_at_qpoint(qpoint=self.dojo_qpt, asr=2, chneut=1, dipdip=1, verbose=1)
        noasr_phbands = ddb.anaget_phmodes_at_qpoint(qpoint=self.dojo_qpt, asr=0, chneut=1, dipdip=1, verbose=1)

        # Convert to JSON and add results to the dojo report.
        # Convert phfreqs[nq, 3* natom] to 1d vector because nq == 1 in this case.
        entry = dict(ecut=self.ecut, pawecutdg=self.dojo_pawecutdg, kppa=self.dojo_kppa,
                     asr2_phfreqs_mev=(asr2_phbands.phfreqs * 1000).ravel().tolist(),
                     noasr_phfreqs_mev=(noasr_phbands.phfreqs * 1000).ravel().tolist()
                )

        self.add_entry_to_dojoreport(entry)
        return results


class RocksaltRelaxationFactory(object):
    """
    Factory producing work objects for the structural relaxation of
    lantanide + nitrogen in rocksalt structure
    """
    def __init__(self, xc):
        """xc is the exchange-correlation functional e.g. PBE, PW."""
        self.xc = XcFunc.asxc(xc)

    def work_for_pseudo(self, pseudo, ecut_list, pawecutdg=None, ngkpt=(8, 8, 8),
                        include_soc=False, workdir=None, manager=None):
        """
        Returns a :class:`Work` object from the given pseudopotential.

        Args:
            pseudo: :class:`Pseudo` object.
            ecut_list: List of cutoff energies in Hartree
            pawecutdg: Cutoff energy of the fine grid (PAW only)
            include_soc: True to include SOC.
            workdir: Working directory.
            manager: :class:`TaskManager` object.
        """
        if pseudo.ispaw and pawecutdg is None:
            raise ValueError("pawecutdg must be specified for PAW calculations.")

        if pseudo.xc != self.xc:
            raise ValueError(
                "Pseudo xc differs from the XC used to instantiate the factory\n"
                "Pseudo: %s, Factory: %s" % (pseudo.xc, self.xc))

        if not ecut_list:
            return None

        # Build rocksalt structure.
        db = raren_database(pseudo.xc)
        a = db.table["ref"][pseudo.symbol]
        species = [pseudo.symbol, "N"]
        structure = Structure.rocksalt(a, species)
        n_pseudo = db.get_n_pseudo(pseudo)

        # Build input template.
        template = abilab.AbinitInput(structure, pseudos=[pseudo, n_pseudo])

        # Input for EuN Rocksalt
        template.set_vars(
            occopt=7,
            tsmear=0.001,
            ecutsm=0.5,
            #kptopt=1,
            ngkpt=ngkpt,
            nshiftk=4,
            shiftk=[0.0, 0.0, 0.5,
                    0.0, 0.5, 0.0,
                    0.5, 0.0, 0.0,
                    0.5, 0.5, 0.5],
            fband=2,
            nstep=60,
            ionmov=2,
            optcell=1,
            dilatmx=1.14,
            tolvrs=1.0E-16,
            tolmxf=1.0E-06,
            ntime=50,
            prtden=0,
            prteig=0,
            prtwf=-1,
        )

        work = RocksaltRelaxationWork()
        work._pseudo = pseudo
        for ecut in ecut_list:
            task = RelaxTask(template.new_with_vars(ecut=ecut))
            #print(task.input)
            work.register(task)

        work.include_soc = include_soc
        return work


class RocksaltRelaxationWork(DojoWork):

    @property
    def dojo_trial(self):
        if not self.include_soc:
            return "raren_relax"
        else:
            return "raren_relax_soc"

    @property
    def dojo_pseudo(self):
        return self._pseudo

    def on_all_ok(self):
        """
        Results are written to the dojoreport.
        """
        def vol2a(vol):
            """Function to compute cubic a0 from primitive v0 (depends on struct_type)"""
            return (4 * vol) ** (1/3.)

        entries = {}
        for task in self:
            ecut = task.input["ecut"]
            #final_structure = task.get_final_structure()
            with task.open_hist() as hist:
                final_structure = hist.final_structure
                initial_energy = hist.etotals[0]

                # Convert float to string with 1 decimal digit.
                dojo_ecut = "%.1f" % ecut
                entries[dojo_ecut] = {
                        "relaxed_a": vol2a(final_structure.volume),
                        "initial_energy_ev_per_atom": float(initial_energy) / len(final_structure),
                }

        #print(entries)
        # Convert to JSON and add results to the dojo report.
        #entry = dict(ecut=self.ecut, pawecutdg=self.dojo_pawecutdg, kppa=self.dojo_kppa)
        #self.add_entry_to_dojoreport(entry)
        #return results
        djrepo = self.djrepo_path

        # Update file content with Filelock.
        with FileLock(djrepo):
            # Read report from file.
            file_report = DojoReport.from_file(djrepo)

            # Create new entry if not already there
            dojo_trial = self.dojo_trial

            #if pop_trial:
            #    file_report.pop(dojo_trial, None)

            if dojo_trial not in file_report:
                file_report[dojo_trial] = {}

            # Convert float to string with 1 decimal digit.
            #dojo_ecut = "%.1f" % self.ecut
            # Check that we are not going to overwrite data.
            #if dojo_ecut in file_report[dojo_trial]:
            #    if not overwrite_data:
            #        raise RuntimeError("dojo_ecut %s already exists in %s. Cannot overwrite data" %
            #                (dojo_ecut, dojo_trial))
            #    else:
            #        file_report[dojo_trial].pop(dojo_ecut)

            # Update file_report by adding the new entry and write new file
            for dojo_ecut, entry in entries.items():
                file_report[dojo_trial][dojo_ecut] = entry

            # Write new dojo report and update the pseudo attribute
            file_report.json_write()
            self._pseudo.dojo_report = file_report

        return dict(returncode=0, message="Lattice paramenters computed and stored in djrepo file")


class RelaxWithGbrvParamsWork(Work):

    def __init__(self, a_guess, struct_type, pseudo, ecut_list=None, pawecutdg=None, ngkpt=(8, 8, 8),
                 spin_mode="unpolarized", include_soc=False, tolvrs=1.e-10, smearing="fermi_dirac:0.001 Ha",
                 ecutsm=0.05, chksymbreak=0, workdir=None, manager=None):
        """
        Build a :class:`Work` for the computation of the relaxed lattice parameter.

        Args:
            structure_type: fcc, bcc
            pseudo: :class:`Pseudo` object.
            ecut_list: Cutoff energy in Hartree
            ngkpt: MP divisions.
            spin_mode: Spin polarization mode.
            toldfe: Tolerance on the energy (Ha)
            smearing: Smearing technique.
            workdir: String specifing the working directory.
            manager: :class:`TaskManager` responsible for the submission of the tasks.
        """
        super(RelaxWithGbrvParamsWork, self).__init__(workdir=workdir, manager=manager)
        self_pseudo = pseudo
        self.include_soc = include_soc
        self.struct_type = struct_type

        if struct_type == "bcc":
            structure = Structure.bcc(a_guess, species=[pseudo.symbol])
        elif struct_type == "fcc":
            structure = Structure.fcc(a_guess, species=[pseudo.symbol])

        # Set extra_abivars.
        extra_abivars = dict(
            pawecutdg=pawecutdg,
            tolvrs=tolvrs,
            prtwf=-1,
            fband=3.0,
            nstep=100,
            ntime=50,
            ecutsm=ecutsm,
            dilatmx=1.1,
        )

        self.ecut_list = ecut_list
        smearing = Smearing.as_smearing(smearing)

        # Kpoint sampling: shiftk depends on struct_type
        shiftk = {"fcc": [0, 0, 0], "bcc": [0.5, 0.5, 0.5]}.get(struct_type)
        spin_mode = SpinMode.as_spinmode(spin_mode)
        ksampling = KSampling.monkhorst(ngkpt, chksymbreak=chksymbreak, shiftk=shiftk,
                                        use_time_reversal=spin_mode.nspinor == 1)
        relax_algo = RelaxationMethod.atoms_and_cell()

        inp = abilab.AbinitInput(structure, pseudo)
        inp.add_abiobjects(ksampling, relax_algo, spin_mode, smearing)
        inp.set_vars(extra_abivars)

        # Register structure relaxation task.
        for ecut in self.ecut_list:
            self.relax_task = self.register_relax_task(inp.new_with_vars(ecut=ecut))

    def on_all_ok(self):
        """
        This method is called when self reaches S_OK. It reads the optimized structure
        from the netcdf file and builds a new work for the computation of the EOS
        with the GBRV parameters.
        """
        # Function to compute cubic a0 from primitive v0 (depends on struct_type)
        vol2a = {"fcc": lambda vol: (4 * vol) ** (1/3.),
                 "bcc": lambda vol: (2 * vol) ** (1/3.),
                 "rocksalt": lambda vol: (4 * vol) ** (1/3.),
                 "ABO3": lambda vol: vol ** (1/3.),
                 "hH": lambda vol: (4 * vol) ** (1/3.),
                 }[self.struct_type]

        results = []
        for task, ecut in zip(self, self.ecut_list):
            structure = task.get_final_structure()
            a0 = vol2a(structure.volume)
            #print("ecut:", ecut, "a0:", a0)
            results.append(dict(ecut=ecut, a0=a0))

        import json
        with open(self.outdir.path_in("a0.json"), "wt") as fh:
            json.dump(results, fh, indent=4) #, sort_keys=True, cls=MontyEncoder)

        return super(RelaxWithGbrvParamsWork, self).on_all_ok()
