# coding: utf-8
"""Base class for Dojo Workkflows."""
from __future__ import division, print_function, unicode_literals

import abc
import sys
import os
import numpy as np
from abipy import abilab

from monty.collections import AttrDict
from monty.pprint import pprint_table
from pymatgen.core.units import Ha_to_eV
from pymatgen.io.abinitio.strategies import ScfStrategy, RelaxStrategy
from pymatgen.io.abinitio.eos import EOS
from pymatgen.io.abinitio.pseudos import Pseudo
from pymatgen.io.abinitio.abiobjects import SpinMode, Smearing, KSampling, Electrons, RelaxationMethod
from pymatgen.io.abinitio.works import Work, build_oneshot_phononwork, OneShotPhononWork
from abipy.core.structure import Structure
from pseudo_dojo.refdata.gbrv import gbrv_database
from pseudo_dojo.refdata.deltafactor import df_database, df_compute

import logging
logger = logging.getLogger(__name__)


class DojoWork(Work):
    __metaclass__ = abc.ABCMeta

    @abc.abstractproperty
    def pseudo(self):
        """:class:`Pseudo` object"""

    @abc.abstractproperty
    def dojo_trial(self):
        """String identifying the DOJO trial. Used to write results in the DOJO_REPORT."""

    def write_dojo_report(self, report, overwrite_data=False):
        """Write/update the DOJO_REPORT section of the pseudopotential."""
        # Read old_report from pseudo.
        old_report = self.pseudo.read_dojo_report()

        dojo_trial = self.dojo_trial
        if dojo_trial not in old_report:
            # Create new entry
            old_report[dojo_trial] = {}
        #else:
        #    # Check that we are not going to overwrite data.
        #    if self.dojo_accuracy in old_report[dojo_trial] and not overwrite_data:
        #        raise RuntimeError("%s already exists in DOJO_REPORT. Cannot overwrite data" % dojo_trial)

        # Update old report card with the new one and write new report
        dojo_ecut = "%.1f" % self.ecut
        old_report[dojo_trial][dojo_ecut] = report
        try:
            self.pseudo.write_dojo_report(old_report)
        except:
            print("Something wrong in write_dojo_report")


def check_conv(values, tol, min_numpts=1, mode="abs", vinf=None):
    """
    Given a list of values and a tolerance tol, returns the leftmost index for which

        abs(value[i] - vinf) < tol if mode == "abs"
    or
        abs(value[i] - vinf) / vinf < tol if mode == "rel"

    returns -1 if convergence is not achieved. By default, vinf = values[-1]

    Args:
        tol: Tolerance
        min_numpts: Minimum number of points that must be converged.
        mode: "abs" for absolute convergence, "rel" for relative convergence.
        vinf: Used to specify an alternative value instead of values[-1].
    """
    vinf = values[-1] if vinf is None else vinf

    if mode == "abs":
        vdiff = [abs(v - vinf) for v in values]
    elif mode == "rel":
        vdiff = [abs(v - vinf) / vinf for v in values]
    else:
        raise ValueError("Wrong mode %s" % mode)

    numpts, i = len(vdiff), -2
    if numpts > min_numpts and vdiff[-2] < tol:
        for i in range(numpts-1, -1, -1):
            if vdiff[i] > tol:
                break
        if (numpts - i -1) < min_numpts: i = -2

    return i + 1


def compute_hints(ecuts, etotals, atols_mev, min_numpts=1, stream=sys.stdout):
    de_low, de_normal, de_high = [a / (1000 * Ha_to_eV) for a in atols_mev]

    num_ene = len(etotals)
    etotal_inf = etotals[-1]

    ihigh = check_conv(etotals, de_high, min_numpts=min_numpts)
    inormal = check_conv(etotals, de_normal)
    ilow = check_conv(etotals, de_low)

    accidx = {"H": ihigh, "N": inormal, "L": ilow}

    table = []; app = table.append

    app(["iter", "ecut", "etotal", "et-e_inf [meV]", "accuracy",])
    for idx, (ec, et) in enumerate(zip(ecuts, etotals)):
        line = "%d %.1f %.7f %.3f" % (idx, ec, et, (et-etotal_inf) * Ha_to_eV * 1.e+3)
        row = line.split() + ["".join(c for c,v in accidx.items() if v == idx)]
        app(row)

    if stream is not None:
        pprint_table(table, out=stream)

    ecut_high, ecut_normal, ecut_low = 3 * (None,)
    exit = (ihigh != -1)

    if exit:
        ecut_low = ecuts[ilow]
        ecut_normal = ecuts[inormal]
        ecut_high = ecuts[ihigh]

    aug_ratios = [1,]
    aug_ratio_low, aug_ratio_normal, aug_ratio_high = 3 * (1,)

    #if not monotonic(etotals, mode="<", atol=1.0e-5):
    #    logger.warning("E(ecut) is not decreasing")
    #    wf_results.push_exceptions("E(ecut) is not decreasing:\n" + str(etotals))

    return AttrDict(
        exit=ihigh != -1,
        etotals=list(etotals),
        ecuts=list(ecuts),
        aug_ratios=aug_ratios,
        low={"ecut": ecut_low, "aug_ratio": aug_ratio_low},
        normal={"ecut": ecut_normal, "aug_ratio": aug_ratio_normal},
        high={"ecut": ecut_high, "aug_ratio": aug_ratio_high})


class PseudoConvergence(DojoWork):

    def __init__(self, pseudo, ecut_slice, nlaunch, atols_mev,
                 toldfe=1.e-8, spin_mode="polarized", acell=(8, 9, 10), 
                 smearing="fermi_dirac:0.1 eV", max_niter=300, workdir=None, manager=None):
        """
        Args:
            pseudo: string or :class:`Pseudo` instance
            ecut_slice: List of cutoff energies or slice object (mainly used for infinite iterations).
            nlaunch:
            atols_mev: List of absolute tolerances in meV (3 entries corresponding to accuracy ["low", "normal", "high"]
            spin_mode: Defined how the electronic spin will be treated.
            acell: Lengths of the periodic box in Bohr.
            smearing: :class:`Smearing` instance or string in the form "mode:tsmear". Default: FemiDirac with T=0.1 eV
            max_niter:
            workdir: Working directory.
            manager: :class:`TaskManager` object.
        """
        super(PseudoConvergence, self).__init__(workdir, manager)

        self._pseudo = Pseudo.as_pseudo(pseudo)
        self.nlaunch = nlaunch; assert nlaunch > 0
        self.atols_mev = atols_mev
        self.toldfe = toldfe
        self.spin_mode = spin_mode
        self.acell = acell
        self.smearing = smearing
        self.max_niter = max_niter; assert max_niter > 0
        self.ecut_slice = ecut_slice; assert isinstance(ecut_slice, slice)

        self.ecuts = []

        if self.pseudo.ispaw:
            raise NotImplementedError("PAW convergence tests are not supported yet")

        for i in range(self.nlaunch):
            ecut = ecut_slice.start + i * ecut_slice.step
            #if self.ecut_slice.stop is not None and ecut > self.ecut_slice.stop: continue
            self.add_task_with_ecut(ecut)

    @property
    def pseudo(self):
        return self._pseudo
                            
    @property
    def dojo_trial(self):
        return "hints"

    def add_task_with_ecut(self, ecut):
        """Register a new task with cutoff energy ecut."""
        # One atom in a box of lenghts acell.
        boxed_atom = Structure.boxed_atom(self.pseudo, acell=self.acell)

        # Gamma-only sampling.
        gamma_only = KSampling.gamma_only()

        extra_abivars = {
            "ecut": ecut,
            "toldfe": self.toldfe,
            "prtwf": 1,
            #"intxc": 1,
        }

        strategy = ScfStrategy(boxed_atom, self.pseudo, gamma_only,
                               spin_mode=self.spin_mode, smearing=self.smearing,
                               **extra_abivars)

        self.ecuts.append(ecut)
        self.register_scf_task(strategy)

    def make_report(self):
        """
        "hints": {
            "high": {"aug_ratio": 1, "ecut": 45},
            "low": {...},
            "normal": {...}
        """
        results = self.work.get_results()
        d = {key: results[key] for key in ["low", "normal", "high"]}

        d.update(dict(
            ecuts=results["ecuts"],
            etotals=results["etotals"]))
        
        if results.exceptions:
            d["_exceptions"] = str(results.exceptions)

        return {self.dojo_key: d}

    def on_all_ok(self):
        """
        This method is called when self reaches S_OK.
        It checks if Etotal(ecut) is converged withing atols_mev
        If the condition is not fulfilled, the callback creates
        nlaunch new tasks with larger values of ecut and we keep on running.
        """
        etotals = self.read_etotals()
        data = compute_hints(self.ecuts, etotals, self.atols_mev)

        if data.exit:
            logger.info("Converged")
            d = {key: data[key] for key in ["low", "normal", "high"]}
                                                                         
            d.update(dict(
                ecuts=data["ecuts"],
                etotals=data["etotals"]))

            #if results.exceptions:
            #    d["_exceptions"] = str(results.exceptions)

            # Read old report from pseudo and add hints
            report = self.pseudo.read_dojo_report()
            report["hints"] = d

            # Write new report
            self.pseudo.write_dojo_report(report)

        else:
            logger.info("Building new tasks")

            estart = self.ecuts[-1] 
            for i in range(self.nlaunch):
                ecut = estart + (i+1) * self.ecut_slice.step
                #if self.ecut_slice.stop is not None and ecut > self.ecut_slice.stop: continue
                self.add_task_with_ecut(ecut)

            if len(self.ecuts) > self.max_niter:
                raise self.Error("Cannot create more that %d tasks, aborting now" % self.max_niter)

            self._finalized = False
            self.flow.allocate()
            self.flow.build_and_pickle_dump()

        return super(PseudoConvergence, self).on_all_ok()


class PPConvergenceFactory(object):
    """
    Factory object that constructs works for analyzing the converge of pseudopotentials.
    """
    def work_for_pseudo(self, pseudo, ecut_slice, nlaunch,
                        toldfe=1.e-8, atols_mev=(10, 1, 0.1), spin_mode="polarized",
                        acell=(8, 9, 10), smearing="fermi_dirac:0.1 eV", workdir=None, manager=None):
        """
        Return a `Work` object given the pseudopotential pseudo.

        Args:
            pseudo: filepath or :class:`Pseudo` object.
            ecut_slice: cutoff energies in Ha units (accepts lists or slice objects)
            toldfe: Tolerance on the total energy (Ha).
            atols_mev: Tolerances in meV for accuracy in ["low", "normal", "high"]
            spin_mode: Spin polarization.
            acell: Length of the real space lattice (Bohr units)
            smearing: Smearing technique.
            workdir: Working directory.
            manager: :class:`TaskManager` object.
        """
        return PseudoConvergence(
            pseudo, ecut_slice, nlaunch, atols_mev,
            toldfe=toldfe, spin_mode=spin_mode,
            acell=acell, smearing=smearing, workdir=workdir, manager=manager)


class DeltaFactoryError(Exception):
    """Base Error class."""


class DeltaFactory(object):
    """Factory class producing work objects for the computation of the delta factor."""
    Error = DeltaFactoryError

    def __init__(self):
        # Get a reference to the deltafactor database
        self._dfdb = df_database()

    def get_cif_path(self, symbol):
        """Returns the path to the CIF file associated to the given symbol."""
        try:
            return self._dfdb.get_cif_path(symbol)
        except KeyError:
            raise self.Error("%s: cannot find CIF file for symbol" % symbol)

    def work_for_pseudo(self, pseudo, accuracy="normal", kppa=6750, ecut=None, pawecutdg=None,
                        toldfe=1.e-9, smearing="fermi_dirac:0.1 eV", workdir=None, manager=None, **kwargs):
        """
        Returns a :class:`Work` object from the given pseudopotential.

        Args:
            kwargs: Extra variables passed to Abinit.

        .. note::

            0.001 Rydberg is the value used with WIEN2K
        """
        pseudo = Pseudo.as_pseudo(pseudo)
        symbol = pseudo.symbol

        if pseudo.ispaw and pawecutdg is None:
            raise ValueError("pawecutdg must be specified for PAW calculations.")

        try:
            cif_path = self.get_cif_path(symbol)
        except Exception as exc:
            raise self.Error(str(exc))

        # Include spin polarization for O, Cr and Mn (antiferromagnetic)
        # and Fe, Co, and Ni (ferromagnetic).
        # antiferromagnetic Cr, O
        # ferrimagnetic Mn
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
            spin_mode = "afm"
            if symbol == 'O':
                kwargs['spinat'] = [(0, 0, 1.5), (0, 0, 1.5), (0, 0, -1.5), (0, 0, -1.5)]
            elif symbol == 'Cr':
                kwargs['spinat'] = [(0, 0, 1.5), (0, 0, -1.5)]
            elif symbol == 'Mn':
                kwargs['spinat'] = [(0, 0, 2.0), (0, 0, 1.9), (0, 0, -2.0), (0, 0, -1.9)]

        # DO NOT CHANGE THE STRUCTURE REPORTED IN THE CIF FILE.
        structure = Structure.from_file(cif_path, primitive=False)

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
            pseudo: String with the name of the pseudopotential file or :class:`Pseudo` object.
            kppa: Number of k-points per atom.
            connect: True if the SCF run should be initialized from the previous run.
            spin_mode: Spin polarization mode.
            toldfe: Tolerance on the energy (Ha)
            smearing: Smearing technique.
            workdir: String specifing the working directory.
            manager: :class:`TaskManager` responsible for the submission of the tasks.
        """
        super(DeltaFactorWork, self).__init__(workdir=workdir, manager=manager)

        self._pseudo = Pseudo.as_pseudo(pseudo)

        spin_mode = SpinMode.as_spinmode(spin_mode)

        # Compute the number of bands from the pseudo and the spin-polarization.
        # Add 6 bands to account for smearing.
        #nval = structure.num_valence_electrons(self.pseudo)
        #spin_fact = 2 if spin_mode.nsppol == 2 else 1
        #nband = int(nval / spin_fact) + 6

        # Set extra_abivars
        self.ecut, self.pawecutdg = ecut, pawecutdg

        extra_abivars = dict(
            ecut=ecut,
            pawecutdg=pawecutdg,
            ecutsm=ecutsm,
            toldfe=toldfe,
            #nband=nband,
            prtwf=0 if not connect else 1,
            #paral_kgb=paral_kgb,
            chkprim=0,
            nstep=200,
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

            ksampling = KSampling.automatic_density(new_structure, kppa, chksymbreak=chksymbreak)

            scf_input = ScfStrategy(new_structure, self.pseudo, ksampling,
                                    accuracy=accuracy, spin_mode=spin_mode,
                                    smearing=smearing, **extra_abivars)

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

    def get_results(self):
        results = super(DeltaFactorWork, self).get_results()

        num_sites = self._input_structure.num_sites
        etotals = self.read_etotals(unit="eV")

        results.update(dict(
            etotals=list(etotals),
            volumes=list(self.volumes),
            num_sites=num_sites))

        d = {}
        try:
            # Use same fit as the one employed for the deltafactor.
            eos_fit = EOS.DeltaFactor().fit(self.volumes/num_sites, etotals/num_sites)

            # Get reference results (Wien2K).
            wien2k = df_database().get_entry(self.pseudo.symbol)

            # Compute deltafactor estimator.
            dfact = df_compute(wien2k.v0, wien2k.b0_GPa, wien2k.b1,
                               eos_fit.v0, eos_fit.b0_GPa, eos_fit.b1, b0_GPa=True)

            dfactprime_meV = dfact * (30 * 100) / (eos_fit.v0 * eos_fit.b0_GPa)

            print("delta", eos_fit)
            print("Ecut %.1f, dfact = %.3f meV, dfactprime %.3f meV" % (self.ecut, dfact, dfactprime_meV))

            results.update({
                "dfact_meV": dfact,
                "v0": eos_fit.v0,
                "b0": eos_fit.b0,
                "b0_GPa": eos_fit.b0_GPa,
                "b1": eos_fit.b1,
                "dfactprime_meV": dfactprime_meV
            })

            d = {k: results[k] for k in 
                ("dfact_meV", "v0", "b0", "b0_GPa", "b1", "etotals", "volumes", "num_sites", "dfactprime_meV")}

            # Write data for the computation of the delta factor
            with open(self.outdir.path_in("deltadata.txt"), "w") as fh:
                fh.write("# Deltafactor = %s meV\n" % dfact)
                fh.write("# Volume/natom [Ang^3] Etotal/natom [eV]\n")
                for v, e in zip(self.volumes, etotals):
                    fh.write("%s %s\n" % (v/num_sites, e/num_sites))

        except EOS.Error as exc:
            results.push_exceptions(exc)

        if results.exceptions:
            d["_exceptions"] = str(results.exceptions)

        self.write_dojo_report(d)

        return results

    def on_all_ok(self):
        """Callback executed when all tasks in self have reached S_OK."""
        return self.get_results()


class GbrvFactory(object):
    """Factory class producing :class:`Work` objects for GBRV calculations."""
    def __init__(self):
        self._db = gbrv_database()

    def make_ref_structure(self, symbol, struct_type, ref):
        """
        Return the structure used in the GBRV tests given the chemical symbol, the structure type
        and the reference code.
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

    def relax_and_eos_work(self, pseudo, struct_type, ecut=None, pawecutdg=None, ref="ae", **kwargs):
        """
        Returns a :class:`Work` object from the given pseudopotential.

        Args:
            kwargs: Extra variables passed to Abinit.

        .. note::

            GBRV tests are done with the following parameteres:

                - No spin polarization for structural relaxation 
                  (only for magnetic moments for which spin-unpolarized structures are used)
                - All calculations are done on an 8x8x8 k-point density and with 0.002 Ry Fermi-Dirac smearing
        """
        pseudo = Pseudo.as_pseudo(pseudo)

        if pseudo.ispaw and pawecutdg is None:
            raise ValueError("pawecutdg must be specified for PAW calculations.")

        structure = self.make_ref_structure(pseudo.symbol, struct_type=struct_type, ref=ref)
 
        return GbrvRelaxAndEosWork(structure, struct_type, pseudo,
                                   ecut=ecut, pawecutdg=pawecutdg, **kwargs)


def gbrv_nband(pseudo):
    # nband/fband are usually too small for the GBRV calculations.
    # FIXME this is not optimal
    nband = pseudo.Z_val
    nband += 0.5 * nband
    nband = int(nband)
    nband = max(nband,  8)
    #print("nband", nband)
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
            pseudo: String with the name of the pseudopotential file or :class:`Pseudo` object.
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
        self._pseudo = Pseudo.as_pseudo(pseudo)
        nband = gbrv_nband(self.pseudo)

        # Set extra_abivars.
        self.extra_abivars = dict(
            ecut=ecut,
            pawecutdg=pawecutdg,
            toldfe=toldfe,
            prtwf=0,
            #ecutsm=0.5,
            nband=nband,
            #paral_kgb=paral_kgb
        )
                                       
        self.extra_abivars.update(**kwargs)
        self.ecut = ecut
        self.smearing = smearing

        # Kpoint sampling: shiftk depends on struct_type
        shiftk = {"fcc": [0, 0, 0], "bcc": [0.5, 0.5, 0.5]}.get(struct_type)
        #ngkpt = (1,1,1)
        self.ksampling = KSampling.monkhorst(ngkpt, chksymbreak=chksymbreak, shiftk=shiftk)
        self.spin_mode = spin_mode
        relax_algo = RelaxationMethod.atoms_and_cell()

        self.relax_input = RelaxStrategy(structure, pseudo, self.ksampling, relax_algo, 
                                         accuracy=accuracy, spin_mode=spin_mode, smearing=smearing, **self.extra_abivars)

        # Register structure relaxation task.
        self.relax_task = self.register_relax_task(self.relax_input)

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
        self.relaxed_structure = relaxed_structure = self.relax_task.read_final_structure()

        # GBRV use nine points from -1% to 1% of the initial guess and fitting the results to a parabola.
        # Note that it's not clear to me if they change the volume or the lattice parameter!
        self.volumes = relaxed_structure.volume * np.arange(99, 101.25, 0.25) / 100.

        for vol in self.volumes:
            new_lattice = relaxed_structure.lattice.scale(vol)
            new_structure = Structure(new_lattice, relaxed_structure.species, relaxed_structure.frac_coords)

            # Add ecutsm
            extra = self.extra_abivars.copy() 
            extra["ecutsm"] = 0.5

            scf_input = ScfStrategy(new_structure, self.pseudo, self.ksampling,
                                    accuracy=self.accuracy, spin_mode=self.spin_mode,
                                    smearing=self.smearing, **extra)

            # Register new task
            self.register_scf_task(scf_input)

        # Allocate new tasks and update the pickle database.
        self.flow.allocate()
        self.flow.build_and_pickle_dump()

    def compute_eos(self):
        results = self.get_results()

        # Read etotals and fit E(V) with a parabola to find minimum
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

        self.write_dojo_report(d)

        return results

    @property
    def add_eos_done(self):
        return len(self) > 1

    def on_all_ok(self):
        """
        This method is called when self reaches S_OK. It reads the optimized structure from the netcdf file and builds
        a new work for the computation of the EOS with the GBRV parameters.
        """
        if not self.add_eos_done:
            logger.info("Building EOS tasks")
            self.add_eos_tasks()
            self._finalized = False
        else:
            logger.info("Computing EOS")
            self.compute_eos()

        return super(GbrvRelaxAndEosWork, self).on_all_ok()


class DFPTError(Exception):
    """Base Error class."""


class DFPTPhononFactory(object):
    """
    Factory class producing `Workflow` objects for DFPT Phonon calculations.
    In particular to test if the acoustic modes are zero
    """

    Error = DFPTError

    def __init__(self, manager=None, workdir=None):
        # reference to the deltafactor database
        # use the elemental solid in the gs configuration
        self._dfdb = df_database()
        self.manager = manager
        self.workdir = workdir

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
        """
        qpoints = [0.00000000E+00,  0.00000000E+00,  0.00000000E+00]
        qpoints = np.reshape(qpoints, (-1, 3))

        # Global variables used both for the GS and the DFPT run.

        kppa = kwargs.pop('kppa')
        ksampling = KSampling.automatic_density(structure, kppa, chksymbreak=0)
        try:
            kwargs.pop('accuracy')
        except:
            pass
        kwargs.pop('smearing')
        global_vars = dict(ksampling.to_abivars(), tsmear=0.005, occopt=7, nstep=200, ecut=12.0, paral_kgb=0)
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
        #global_vars.pop('#comment')
        electrons = structure.num_valence_electrons(pseudos)
        global_vars.update(nband=electrons)
        global_vars.update(nbdbuf=int(electrons/4))

        inp = abilab.AbiInput(pseudos=pseudos, ndtset=1+len(qpoints))
        inp.set_structure(structure)
        inp.set_variables(**global_vars)

        for i, qpt in enumerate(qpoints):
            # Response-function calculation for phonons.
            inp[i+2].set_variables(nstep=200, iscf=7, rfphon=1, nqpt=1, qpt=qpt, kptopt=2, rfasr=2)

        # Split input into gs_inp and ph_inputs
        return inp.split_datasets()

    def work_for_pseudo(self, pseudo, **kwargs):
        """
        Create an `AbinitFlow` for phonon calculations:

            1) One workflow for the GS run.

            2) nqpt workflows for phonon calculations. Each workflow contains
               nirred tasks where nirred is the number of irreducible phonon perturbations
               for that particular q-point.
        """
        accuracy = kwargs.pop('accuracy')

        pseudo = Pseudo.as_pseudo(pseudo)

        structure_or_cif = self.get_cif_path(pseudo.symbol)

        if not isinstance(structure_or_cif, Structure):
            # Assume CIF file
            structure = Structure.from_file(structure_or_cif, primitive=False)
        else:
            structure = structure_or_cif
        #print(structure)

        all_inps = self.scf_ph_inputs(pseudos=[pseudo], structure=structure, **kwargs)
        scf_input, ph_inputs = all_inps[0], all_inps[1:]

        # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
        workdir = self.workdir
        if not self.workdir:
            workdir = os.path.basename(__file__).replace(".py", "").replace("run_", "flow_")

        # Instantiate the TaskManager.
        manager = abilab.TaskManager.from_user_config() if not self.manager else \
            abilab.TaskManager.from_file(self.manager)

        work = build_oneshot_phononwork(scf_input=scf_input, ph_inputs=ph_inputs, work_class=PhononDojoWorkflow)
        work._pseudo = pseudo
#        work.set_dojo_accuracy(accuracy=accuracy)
        return work


class PhononDojoWorkflow(OneShotPhononWork, DojoWork):
    @property
    def dojo_trial(self):
        return "phonon"

    @property
    def pseudo(self):
        return self._pseudo

#        return super(GbrvRelaxAndEosWork, self).on_all_ok()
