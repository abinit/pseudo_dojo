# coding: utf-8
"""Base class for Dojo Workkflows."""
from __future__ import division, print_function, unicode_literals

import abc
import sys
import os
import json
import logging
import numpy as np

from monty.collections import AttrDict
from monty.pprint import pprint_table
from monty.io import FileLock
from pymatgen.core.xcfunc import XcFunc
from pymatgen.analysis.eos import EOS
from pymatgen.io.abinit.abiobjects import SpinMode, Smearing, KSampling, RelaxationMethod
from pymatgen.io.abinit.works import Work, build_oneshot_phononwork, OneShotPhononWork
from abipy.core.structure import Structure
from abipy import abilab
from pseudo_dojo.core.dojoreport import DojoReport, compute_dfact_entry
from pseudo_dojo.refdata.gbrv import gbrv_database
from pseudo_dojo.refdata.deltafactor import df_database, df_compute


logger = logging.getLogger(__name__)


class DojoWork(Work):
    __metaclass__ = abc.ABCMeta

    @abc.abstractproperty
    def pseudo(self):
        """:class:`Pseudo` object"""

    @abc.abstractproperty
    def dojo_trial(self):
        """String identifying the DOJO trial. Used to write results in the DOJO_REPORT."""

    def add_entry_to_dojoreport(self, entry, overwrite_data=False):
        """
        Write/update the DOJO_REPORT section of the pseudopotential.
        Important paramenters such as the name of the dojo_trial and the energy cutoff
        are provided by the sub-class.
        Client code is responsible for preparing the dictionary with the data.

        Args:
            entry: Dictionary with results.
            overwrite_data: If False, the routine raises an exception if this entry is
                already filled.
        """
        root, ext = os.path.splitext(self.pseudo.filepath)
        djrepo = root + ".djrepo"

        # Update file content with Filelock.
        with FileLock(djrepo):
            # Read report from file.
            file_report = DojoReport.from_file(djrepo)

            # Create new entry if not already there
            dojo_trial = self.dojo_trial
            if dojo_trial not in file_report: file_report[dojo_trial] = {}

            # Convert float to string with 1 decimal digit.
            dojo_ecut = "%.1f" % self.ecut

            # Check that we are not going to overwrite data.
            if dojo_ecut in file_report[dojo_trial] and not overwrite_data:
                raise RuntimeError("dojo_ecut %s already exists in %s. Cannot overwrite data" % (dojo_ecut, dojo_trial))

            # Update file_report by adding the new entry and write new file
            file_report[dojo_trial][dojo_ecut] = entry

            # Write new dojo report and update the pseudo attribute
            file_report.json_write(djrepo)
            self._pseudo.dojo_report = file_report


class FactoryError(Exception):
    """Base Error class raised by Factory objects."""


class EbandsFactory(object):
    """Factroy class producing work objects for the calculation of ebands for testing for ghosts."""

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

    def work_for_pseudo(self, pseudo, accuracy="normal", kppa=3000, ecut=None, pawecutdg=None, spin_mode="unpolarized",
                        toldfe=1.e-9, smearing="fermi_dirac:0.1 eV", workdir=None, manager=None, **kwargs):
        """
        Returns a :class:`Work` object from the given pseudopotential.

        Args:
            pseudo: :class:`Pseudo` object.
            kppa: Number of k-points per atom.
            ecut: Cutoff energy in Hartree
            pawecutdg: Cutoff energy of the fine grid (PAW only)
            spin_mode: Spin polarization option
            toldfe: Tolerance on the total energy (Ha).
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

        return EbandsFactorWork(
            structure, pseudo, kppa,
            spin_mode=spin_mode, toldfe=toldfe, smearing=smearing,
            accuracy=accuracy, ecut=ecut, pawecutdg=pawecutdg, ecutsm=0.5,
            workdir=workdir, manager=manager, **kwargs)


class EbandsFactorWork(DojoWork):
    """Work for the calculation of the deltafactor."""

    # THIS IS WRONG! toldfe for band structure, lot of bands computed with SCF!
    def __init__(self, structure, pseudo, kppa,
                 bands_factor=10, ecut=None, pawecutdg=None, ecutsm=0.5,
                 spin_mode="polarized", toldfe=1.e-9, smearing="fermi_dirac:0.1 eV",
                 accuracy="normal", chksymbreak=0, workdir=None, manager=None, **kwargs):
        """
        Build a :class:`Work` for the computation of a bandstructure to check for ghosts.

        Args:
            structure: :class:`Structure` object
            pseudo: String with the name of the pseudopotential file or :class:`Pseudo` object.
            kppa: Number of k-points per atom.
            bands_factor: Number of bands computed is given by bands_factor * int(nval / nsppol)
            spin_mode: Spin polarization mode.
            toldfe: Tolerance on the energy (Ha)
            smearing: Smearing technique.
            workdir: String specifing the working directory.
            manager: :class:`TaskManager` responsible for the submission of the tasks.
        """
        super(EbandsFactorWork, self).__init__(workdir=workdir, manager=manager)
        self._pseudo = pseudo

        spin_mode = SpinMode.as_spinmode(spin_mode)
        smearing = Smearing.as_smearing(smearing)

        # Compute the number of bands from the pseudo and the spin-polarization.
        # Take 10 times the number of bands to sample the empty space.
        nval = structure.num_valence_electrons(self.pseudo)
        nband = bands_factor * int(nval / spin_mode.nsppol)

        # Set extra_abivars
        self.ecut, self.pawecutdg = ecut, pawecutdg

        extra_abivars = dict(
            ecut=ecut,
            pawecutdg=pawecutdg,
            ecutsm=ecutsm,
            nband=nband,
            toldfe=toldfe,
            prtwf=-1,
            nstep=200,
            chkprim=0,
            mem_test=0
        )

        extra_abivars.update(**kwargs)

        # Disable time-reversal if nspinor == 2
        ksampling = KSampling.automatic_density(structure, kppa, chksymbreak=chksymbreak,
                                                use_time_reversal=spin_mode.nspinor==1)

        scf_input = abilab.AbinitInput(structure=structure, pseudos=self.pseudo)
        scf_input.add_abiobjects(ksampling, smearing, spin_mode)
        scf_input.set_vars(extra_abivars)
        self.register_scf_task(scf_input)

    @property
    def pseudo(self):
        return self._pseudo

    @property
    def dojo_trial(self):
        return "ebands"

    def on_all_ok(self):
        """Callback executed when all tasks in the work have reached S_OK."""
        # store the path of the GSR nc file in the dojoreport
        # during the validation the bandstrure is plotted and the validator is asked to give the energy up to which
        # no sign of ghosts is visible
        # during plot the band structuur is plotted as long as a filename is present and the file can be found if the
        # if the energy is present the energy is given

        #TODO fix magic
        path = str(self.workdir)
        outfile = os.path.join(str(self[0].outdir), "out_GSR.nc")
        entry = {'workdir': path, 'GSR-nc': outfile}

        self.add_entry_to_dojoreport(entry)

        return entry
        #from abipy.abilab import abiopen
        #with abiopen(path) as gsr:
        #    ebands = gsr.ebands
        #    fig = ebands.plot_with_edos(ebands.get_edos(width=0.05, step=0.02))
        #    return fig


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

    def work_for_pseudo(self, pseudo, accuracy="normal", kppa=6750, ecut=None, pawecutdg=None,
                        toldfe=1.e-9, smearing="fermi_dirac:0.1 eV", include_soc=False,
                        workdir=None, manager=None, **kwargs):
        """
        Returns a :class:`Work` object from the given pseudopotential.

        Args:
            pseudo: :class:`Pseudo` object.
            kppa: kpoint per atom
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

        # Include spin polarization for O, Cr and Mn (antiferromagnetic)
        # and Fe, Co, and Ni (ferromagnetic).
        # antiferromagnetic Cr, O, ferrimagnetic Mn
        spin_mode = "unpolarized"

        if symbol in ["Fe", "Co", "Ni"]:
            spin_mode = "polarized"
            if symbol == "Fe":
                kwargs['spinat'] = 2 * [(0, 0, 2.3)]
            if symbol == "Co":
                kwargs['spinat'] = 2 * [(0, 0, 1.2)]
            if symbol == "Ni":
                kwargs['spinat'] = 4 * [(0, 0, 0.6)]

        if symbol in ["O", "Cr", "Mn"]:
            # Here we could have problems with include_so since we don't enforce "afm"
            spin_mode = "afm"
            if symbol == 'O':
                kwargs['spinat'] = [(0, 0, 1.5), (0, 0, 1.5), (0, 0, -1.5), (0, 0, -1.5)]
            elif symbol == 'Cr':
                kwargs['spinat'] = [(0, 0, 1.5), (0, 0, -1.5)]
            elif symbol == 'Mn':
                kwargs['spinat'] = [(0, 0, 2.0), (0, 0, 1.9), (0, 0, -2.0), (0, 0, -1.9)]

        if include_soc: spin_mode = "spinor"

        # Magnetic elements:
        # Start from previous SCF run to avoid getting trapped in local minima
        connect = symbol in ("Fe", "Co", "Ni", "Cr", "Mn", "O", "Zn", "Cu")

        return DeltaFactorWork(
            structure, pseudo, kppa, connect,
            spin_mode=spin_mode, toldfe=toldfe, smearing=smearing,
            accuracy=accuracy, ecut=ecut, pawecutdg=pawecutdg, ecutsm=0.5,
            workdir=workdir, manager=manager, **kwargs)


class DeltaFactorWork(DojoWork):
    """Work for the calculation of the deltafactor."""

    def __init__(self, structure, pseudo, kppa, connect,
                 ecut=None, pawecutdg=None, ecutsm=0.5,
                 spin_mode="polarized", toldfe=1.e-9, smearing="fermi_dirac:0.1 eV",
                 accuracy="normal", chksymbreak=0, workdir=None, manager=None, **kwargs):
        """
        Build a :class:`Work` for the computation of the deltafactor.

        Args:
            structure: :class:`Structure` object
            pseudo: :class:`Pseudo` object.
            kppa: Number of k-points per atom.
            connect: True if the SCF run should be initialized from the previous run.
            spin_mode: Spin polarization mode.
            toldfe: Tolerance on the energy (Ha)
            smearing: Smearing technique.
            workdir: String specifing the working directory.
            manager: :class:`TaskManager` responsible for the submission of the tasks.
        """
        super(DeltaFactorWork, self).__init__(workdir=workdir, manager=manager)
        self._pseudo = pseudo

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
            new_lattice = structure.lattice.scale(vol)

            new_structure = Structure(new_lattice, structure.species, structure.frac_coords)

            ksampling = KSampling.automatic_density(new_structure, kppa, chksymbreak=chksymbreak,
                                                    use_time_reversal=spin_mode.nspinor==1)

            scf_input = abilab.AbinitInput(structure=new_structure, pseudos=self.pseudo)
            scf_input.add_abiobjects(ksampling, smearing, spin_mode)
            scf_input.set_vars(extra_abivars)

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
    def pseudo(self):
        return self._pseudo

    @property
    def dojo_trial(self):
        return "deltafactor"

    def on_all_ok(self):
        """
        This method is called once the `Work` is completed i.e. when all the tasks
        have reached status S_OK. Here we gather the results of the different tasks,
        the deltafactor is computed and the results are stored in the JSON file.
        """
        num_sites = self._input_structure.num_sites
        etotals = self.read_etotals(unit="eV")

        d, eos_fit = compute_dfact_entry(self.pseudo, num_sites, self.volumes, etotals)

        print("[%s]" % self.pseudo.symbol, "eos_fit:", eos_fit)
        print("Ecut %.1f, dfact = %.3f meV, dfactprime %.3f meV" % (self.ecut, d["dfact_meV"], d["dfactprime_meV"]))

        self.add_entry_to_dojoreport(d)

        return dict(returncode=0, message="Delta factor computed")


class GbrvFactory(object):
    """Factory class producing :class:`Work` objects for GBRV calculations."""

    def __init__(self, xc):
        """xc: exchange-correlation functional e.g. PBE or PW."""
        self.xc = XcFunc.asxc(xc)
        if self.xc != "PBE":
            raise ValueError("Gbrv database supports only PBE pseudos")
        self._db = gbrv_database()

    def make_ref_structure(self, symbol, struct_type, ref):
        """
        Return the structure used in the GBRV tests given the chemical symbol,
        the structure type and the reference code.
        """
        # Get the entry in the database
        entry = self._db.get_entry(symbol, struct_type)

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
                                   ecut=ecut, pawecutdg=pawecutdg, spin_mode=spin_mode, **kwargs)


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


class GbrvRelaxAndEosWork(DojoWork):

    def __init__(self, structure, struct_type, pseudo, ecut=None, pawecutdg=None, ngkpt=(8, 8, 8),
                 spin_mode="unpolarized", toldfe=1.e-9, smearing="fermi_dirac:0.001 Ha",
                 accuracy="normal", ecutsm=0.05, chksymbreak=0,
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
        self.accuracy = accuracy

        # nband must be large enough to accomodate fractional occupancies.
        fband = kwargs.pop("fband", None)
        self._pseudo = pseudo
        nband = gbrv_nband(self.pseudo)

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
                                            use_time_reversal=self.spin_mode.nspinor==1)
        relax_algo = RelaxationMethod.atoms_and_cell()

        inp = abilab.AbinitInput(structure, pseudo)
        inp.add_abiobjects(self.ksampling, relax_algo, self.spin_mode, self.smearing)
        inp.set_vars(self.extra_abivars)

        # Register structure relaxation task.
        self.relax_task = self.register_relax_task(inp)

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

            scf_input = abilab.AbinitInput(new_structure, self.pseudo)
            scf_input.add_abiobjects(self.ksampling, self.spin_mode, self.smearing)
            scf_input.set_vars(extra)

            # Register new task
            self.register_scf_task(scf_input)

        # Allocate new tasks and update the pickle database.
        self.flow.allocate()
        self.flow.build_and_pickle_dump()

    def compute_eos(self):
        self.history.info("Computing EOS")

        results = self.get_results()

        # Read etotals and fit E(V) with a parabola to find the minimum
        etotals = self.read_etotals(unit="eV")[1:]
        assert len(etotals) == len(self.volumes)

        results.update(dict(
            etotals=list(etotals),
            volumes=list(self.volumes),
            num_sites=len(self.relaxed_structure),
        ))

        try:
            eos_fit = EOS.Quadratic().fit(self.volumes, etotals)
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
            struct_type=self.struct_type
        ))

        db = gbrv_database()
        entry = db.get_entry(self.pseudo.symbol, stype=self.struct_type)

        pawabs_err = a0 - entry.gbrv_paw
        pawrel_err = 100 * (a0 - entry.gbrv_paw) / entry.gbrv_paw

        # AE results for P and Hg are missing.
        if entry.ae is not None:
            abs_err = a0 - entry.ae
            rel_err = 100 * (a0 - entry.ae) / entry.ae
        else:
            # Use GBRV_PAW as reference.
            abs_err = pawabs_err
            rel_err = pawrel_err

        print("for GBRV struct_type: ", self.struct_type, "a0= ", a0, "Angstrom")
        print("AE - THIS: abs_err = %f, rel_err = %f %%" % (abs_err, rel_err))
        print("GBRV-PAW - THIS: abs_err = %f, rel_err = %f %%" % (pawabs_err, pawrel_err))

        d = {k: results[k] for k in ("a0", "etotals", "volumes")}
        d["a0_abs_err"] = abs_err
        d["a0_rel_err"] = rel_err
        if results.exceptions:
            d["_exceptions"] = str(results.exceptions)

        self.add_entry_to_dojoreport(d)

        return results

    @property
    def add_eos_done(self):
        return len(self) > 1

    def on_all_ok(self):
        """
        This method is called when self reaches S_OK. It reads the optimized structure
        from the netcdf file and builds a new work for the computation of the EOS
        with the GBRV parameters.
        """
        if not self.add_eos_done:
            self.add_eos_tasks()
            self._finalized = False
        else:
            self.compute_eos()

        return super(GbrvRelaxAndEosWork, self).on_all_ok()


class DFPTPhononFactory(object):
    """
    Factory class producing `Workflow` objects for DFPT Phonon calculations.
    The work tests if the acoustic modes at Gamma are zero, or at least from which cuttoff
    they can be made zero by imposing the accoustic sum rule.
    """

    Error = FactoryError

    def __init__(self, xc):
        """xc is the exchange-correlation functional e.g. PBE, PW."""
        # Get reference to the deltafactor database
        # Use the elemental solid in the gs configuration
        self._dfdb = df_database(xc)

    def get_cif_path(self, symbol):
        """Returns the path to the CIF file associated to the given symbol."""
        try:
            return self._dfdb.get_cif_path(symbol)
        except KeyError:
            raise self.Error("%s: cannot find CIF file for symbol" % symbol)

    @staticmethod
    def scf_ph_inputs(structure, pseudos, **kwargs):
        """
        This function constructs the input files for the phonon calculation:
        GS input + the input files for the phonon calculation.

        kwargs:
            ecut: the ecut at which the input is generated
            kppa: kpoint per atom
            smearing: is removed
            qpt: optional, list of qpoints. if not present gamma is added

        the rest are passed as abinit input variables
        """
        qpoints = kwargs.pop('qpt', [0.00000000E+00,  0.00000000E+00,  0.00000000E+00])
        qpoints = np.reshape(qpoints, (-1, 3))

        # Global variables used both for the GS and the DFPT run.
        kppa = kwargs.pop('kppa')
        ksampling = KSampling.automatic_density(structure, kppa, chksymbreak=0)
        try:
            kwargs.pop('accuracy')
        except KeyError:
            pass

        kwargs.pop('smearing')
        # to be applicable to all systems we treat all as is they were metals
        # some systems have a non primitive cell to allow for a anti ferromagnetic structure > chkprim = 0
        global_vars = dict(ksampling.to_abivars(), tsmear=0.005, occopt=7, nstep=200, ecut=12.0, paral_kgb=0, chkprim=0)
        global_vars.update(**kwargs)

        # if not tolwfr is specified explicitly we remove any other tol and put tolwfr = 1e-16
        tolwfr = 1e-20
        for k in global_vars.keys():
            if 'tol' in k:
                if k == 'tolwfr':
                    tolwfr = global_vars.pop(k)
                else:
                    global_vars.pop(k)

        global_vars['tolwfr'] = tolwfr
        electrons = structure.num_valence_electrons(pseudos)
        global_vars.update(nband=electrons)
        global_vars.update(nbdbuf=int(electrons/4))

        multi = abilab.MultiDataset(structure=structure, pseudos=pseudos, ndtset=1+len(qpoints))
        multi.set_vars(global_vars)

        # Check xc:
        #if any(p.xc != self._dbdb.xc for p in multi[0].pseudos):
        #raise ValueError("XC found in pseudos do not agree with the one used in the factory")

        rfasr = kwargs.pop('rfasr', 2)

        for i, qpt in enumerate(qpoints):
            # Response-function calculation for phonons.
            # rfatpol=[1, natom],  # Set of atoms to displace.
            # rfdir=[1, 1, 1],     # Along this set of reduced coordinate axis
            multi[i+1].set_vars(nstep=200, iscf=7, rfphon=1, nqpt=1, qpt=qpt, kptopt=2, rfasr=rfasr,
                                rfatpol=[1, len(structure)], rfdir=[1, 1, 1])

            # rfasr = 1 is not correct
            # response calculations cannot be restarted > nstep = 200, a problem to solve here is that abinit continues
            # happily even is NaN are produced ... TODO fix abinit

        # Split input into gs_inp and ph_inputs
        return multi.split_datasets()

    def work_for_pseudo(self, pseudo, workdir=None, manager=None, **kwargs):
        """
        Create a :class:`Work` for phonon calculations:

            1) One workflow for the GS run.

            2) nqpt workflows for phonon calculations. Each workflow contains
               nirred tasks where nirred is the number of irreducible phonon perturbations
               for that particular q-point.

        Args:
            pseudo: filepath or :class:`Pseudo` object.
            workdir: Working directory.
            manager: :class:`TaskManager` object.

            the kwargs are passed to scf_hp_inputs
        """
        try:
            qpt = kwargs['qpt']
        except IndexError:
            raise ValueError('A phonon test needs to specify a qpoint.')

        kwargs.pop('accuracy')

        if pseudo.xc != self._dfdb.xc:
            raise ValueError(
                "Pseudo xc differs from the XC used to instantiate the factory\n"
                "Pseudo: %s, Database: %s" % (pseudo.xc, self._dfdb.xc))

        structure_or_cif = self.get_cif_path(pseudo.symbol)

        if not isinstance(structure_or_cif, Structure):
            # Assume CIF file
            structure = Structure.from_file(structure_or_cif, primitive=False)
        else:
            structure = structure_or_cif

        nat = len(structure)
        report = pseudo.dojo_report
        ecut_str = '%.1f' % kwargs['ecut']

        try:
            v0 = nat * report['deltafactor'][ecut_str]['v0']
        except KeyError:
            try:
                v0 = nat * report['deltafactor'][float(ecut_str)]['v0']
            except KeyError:
                logger.critical("The deltafactor calculation for this ecut is not done. The phonon task can not be created")
                logger.critical("Returning None")
                return None

        structure.scale_lattice(v0)

        if 'rfasr' in kwargs:
            trial_name = 'phwoa'
        else:
            trial_name = 'phonon'

        all_inps = self.scf_ph_inputs(pseudos=[pseudo], structure=structure, **kwargs)
        scf_input, ph_inputs = all_inps[0], all_inps[1:]

        work = build_oneshot_phononwork(scf_input=scf_input, ph_inputs=ph_inputs, work_class=PhononDojoWork,
                                        workdir=workdir, manager=manager)
        #print('after build_oneshot_phonon', work)
        work.set_dojo_trial(qpt, trial_name)
        #print(scf_input.keys())
        work.ecut = scf_input['ecut']
        work._pseudo = pseudo

        return work


class PhononDojoWork(OneShotPhononWork, DojoWork):
    @property
    def dojo_trial(self):
        return self._trial

    def set_dojo_trial(self, qpt, trial_name):
        if max(qpt) == 0:
            self._trial = trial_name
        elif max(qpt) == 0.5 and min(qpt) == 0.5:
            self._trial = trial_name + '_hhh'
        else:
            raise ValueError('Only dojo phonon works of the Gamma and 0.5, 0.5, 0.5 qpoints have been implemented')

    @property
    def pseudo(self):
        return self._pseudo

    def on_all_ok(self):
        d = self.get_results()
        entry = d['phonons'][0].freqs.tolist()
        self.add_entry_to_dojoreport(entry)
        return d
