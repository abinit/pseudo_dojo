# coding: utf-8
"""Workflows to perform GBRV tests for binary and ternary compunds"""
import numpy as np
import pandas as pd
import json

from abipy import abilab
from abipy.flowtk.abiobjects import SpinMode, Smearing, KSampling, RelaxationMethod
from abipy.core.structure import Structure
from abipy.flowtk.works import Work
from abipy.flowtk.flows import Flow
from pseudo_dojo.core.pseudos import DojoTable
from pseudo_dojo.refdata.gbrv import gbrv_database
from pseudo_dojo.dojo.gbrv_outdb import GbrvOutdb
from pseudo_dojo.util.dojo_eos import EOS

import logging
logger = logging.getLogger(__name__)


class GbrvCompoundsFactory(object):
    """Factory class producing :class:`Work` objects for GBRV calculations."""
    def __init__(self, xc):
        """xc: Exchange-correlation functional."""
        self.db = gbrv_database(xc)

    def make_ref_structure(self, formula, struct_type, ref="ae"):
        """
        Return the structure used in the GBRV tests given the chemical formula,
        the structure type and the reference code.
        """
        # Get the entry in the database
        entry = self.db.get_entry(formula, struct_type)
        if entry is None:
            logger.critical("Cannot find entry for %s, returning None!" % formula)
            return None

        # Build the structure and handle a possibly missing value.
        structure = entry.build_structure(ref=ref)

        if structure is None:
            logger.warning("No AE structure for %s. Will use gbrv_uspp data." % formula)
            structure = entry.build_structure(ref="gbrv_uspp")
        if structure is None:
            logger.critical("Cannot initialize structure for %s, returning None!" % formula)

        return structure

    def relax_and_eos_work(self, accuracy, pseudos, formula, struct_type,
                           ecut=None, pawecutdg=None, ref="ae", ngkpt=(8, 8, 8), fband=2.0):
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
        pseudos = DojoTable.as_table(pseudos)
        if pseudos.allpaw and pawecutdg is None:
            raise ValueError("pawecutdg must be specified for PAW calculations.")

        structure = self.make_ref_structure(formula, struct_type=struct_type, ref=ref)

        return GbrvCompoundRelaxAndEosWork(structure, formula, struct_type, pseudos, self.db.xc, accuracy,
                                           ecut=ecut, pawecutdg=pawecutdg, ngkpt=ngkpt, fband=fband)


class GbrvCompoundRelaxAndEosWork(Work):

    def __init__(self, structure, formula, struct_type, pseudos, xc, accuracy,
                 ecut=None, pawecutdg=None, ngkpt=(8, 8, 8), fband=2.0,
                 spin_mode="unpolarized", smearing="fermi_dirac:0.001 Ha",
                 chksymbreak=0, workdir=None, manager=None):
        """
        Build a :class:`Work` for the computation of the relaxed lattice parameter.

        Args:
            structure: :class:`Structure` object.
            struct_type: fcc, bcc
            pseudos: Pseudopotentials
            xc: Exchange-correlation type.
            accuracy:
            ecut: Cutoff energy in Hartree
            ngkpt: Divisions for k-mesh.
            fband: Input variable, used to compute the number of bands.
            spin_mode: Spin polarization mode.
            smearing: Smearing technique.
            chksymbreak: Input variable.
            workdir: String specifing the working directory.
            manager: :class:`TaskManager` responsible for the submission of the tasks.
        """
        super(GbrvCompoundRelaxAndEosWork, self).__init__(workdir=workdir, manager=manager)

        self.pseudos = pseudos
        self.xc = xc
        if (any(xc != p.xc for p in pseudos)):
            raise ValueError("Input XC does not agree with XC from pseudos.")

        self.formula = formula
        self.struct_type = struct_type
        self.accuracy = accuracy
        self.ecut, self.pawecutdg = ecut, pawecutdg

        # Set extra_abivars.
        # Use tolvrs for the relaxation, toldfe for the EOS
        self.extra_abivars = dict(
            ecut=ecut,
            pawecutdg=pawecutdg,
            tolvrs=1e-10,
            ecutsm=0.5,
            mem_test=0,
            fband=fband,
            # nband must be large enough to accomodate fractional occupancies.
            #paral_kgb=kwargs.pop("paral_kgb", 0),
            #nband=nband,
        )

        self.smearing = Smearing.as_smearing(smearing)

        # Kpoint sampling: shiftk depends on struct_type
        #shiftk = {"fcc": [0, 0, 0], "bcc": [0.5, 0.5, 0.5]}.get(struct_type)
        shiftk = [0, 0, 0]
        self.ksampling = KSampling.monkhorst(ngkpt, chksymbreak=chksymbreak, shiftk=shiftk)
        self.spin_mode = SpinMode.as_spinmode(spin_mode)
        relax_algo = RelaxationMethod.atoms_and_cell()

        inp = abilab.AbinitInput(structure, pseudos)
        inp.add_abiobjects(self.ksampling, relax_algo, self.spin_mode, self.smearing)
        inp.set_vars(self.extra_abivars)

        # Register structure relaxation task.
        self.relax_task = self.register_relax_task(inp)

    def add_eos_tasks(self):
        """
        Read the optimized structure from the netcdf file and add to self a new
        a new list of `ScfTask` for the computation of the EOS with the GBRV setup.
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

            scf_input = abilab.AbinitInput(new_structure, self.pseudos)

            scf_input.set_vars(self.extra_abivars.copy())
            scf_input.add_abiobjects(self.ksampling, self.spin_mode, self.smearing)
            # Use toldfe instead of tolvrs
            scf_input.pop_tolerances()
            scf_input.set_vars(toldfe=1e-10)

            # Register new task
            self.register_scf_task(scf_input)

        # Allocate new tasks and update the pickle database.
        self.flow.allocate()
        self.flow.build_and_pickle_dump()

    def compute_eos(self):
        """Compute the EOS as described in the GBRV paper."""
        self.history.info("Computing EOS")

        # Read etotals and fit E(V) with a parabola to find the minimum
        etotals = self.read_etotals(unit="eV")[1:]
        assert len(etotals) == len(self.volumes)

        results = {}
        results.update(dict(
            etotals=list(etotals),
            volumes=list(self.volumes),
            num_sites=len(self.relaxed_structure),
            #pseudos=self.pseudos.as_dict(),
        ))

        try:
            eos_fit = EOS.Quadratic().fit(self.volumes, etotals)
        except EOS.Error as exc:
            results.push_exceptions(exc)

        # Function to compute cubic a0 from primitive v0 (depends on struct_type)
        vol2a = {"fcc": lambda vol: (4 * vol) ** (1/3.),
                 "bcc": lambda vol: (2 * vol) ** (1/3.),
                 "rocksalt": lambda vol: (4 * vol) ** (1/3.),
                 "ABO3": lambda vol: vol ** (1/3.),
                 "hH": lambda vol: (4 * vol) ** (1/3.),
                 }[self.struct_type]
        a0 = vol2a(eos_fit.v0)

        results.update(dict(
            ecut=self.ecut,
            pawecutdg=self.pawecutdg,
            struct_type=self.struct_type,
            v0=eos_fit.v0,
            b0=eos_fit.b0,
            #b1=eos_fit.b1, # infinity
            a0=a0,
        ))

        db = gbrv_database(xc=self.xc)
        entry = db.get_entry(self.formula, stype=self.struct_type)

        try:
            pawabs_err = a0 - entry.gbrv_paw
            pawrel_err = 100 * (a0 - entry.gbrv_paw) / entry.gbrv_paw
        except Exception:
            # Some paw_abinit entries are not available.
            pawabs_err = np.inf
            pawrel_err = np.inf

        # If AE results are missing we use GBRV_PAW as reference.
        if entry.ae is not None:
            ref_a0 = entry.ae
            abs_err = a0 - entry.ae
            rel_err = 100 * (a0 - entry.ae) / entry.ae
        else:
            assert pawabs_err is not None
            assert pawrel_err is not None
            ref_a0 = entry.gbrv_paw
            abs_err = pawabs_err
            rel_err = pawrel_err

        results["a0_rel_err"] = rel_err
        results["a0_abs_err"] = abs_err
        results["pseudos"] = {p.basename: p.md5 for p in self.pseudos}

        print(80 * "=")
        print("pseudos:", list(p.basename for p in self.pseudos))
        print("for %s (%s): my_a0=%.3f, ref_a0=%.3f Angstrom" % (self.formula, self.struct_type, a0, ref_a0))
        print("AE - THIS: abs_err=%f, rel_err=%f %%" % (abs_err, rel_err))
        print("GBRV-PAW - THIS: abs_err=%f, rel_err=%f %%" % (pawabs_err, pawrel_err))
        print(80 * "=")

        # Write results in outdir in JSON format.
        with open(self.outdir.path_in("gbrv_results.json"), "wt") as fh:
            json.dump(results, fh, indent=-1, sort_keys=True)

        # Update the database.
        if self.outdb_path is not None:
            GbrvOutdb.insert_results(self.outdb_path, self.struct_type, self.formula,
                                     self.accuracy, self.pseudos, results)

        return results

    @property
    def add_eos_done(self):
        return len(self) > 1

    def on_all_ok(self):
        """
        This method is called when self reaches S_OK. It reads the optimized
        structure from the netcdf file and builds a new work for the computation
        of the EOS with the GBRV parameters.
        """
        if not self.add_eos_done:
            self.add_eos_tasks()
            self.finalized = False
        else:
            self.compute_eos()

        return super(GbrvCompoundRelaxAndEosWork, self).on_all_ok()

    def set_outdb(self, path):
        """
        This function set the outdb property (a database with the Gbrv results for the compounds)
        Use this function when you want the work to write the results of the
        calculation to the ouddb calculation.
        """
        self._outdb_path = path
        return self

    @property
    def outdb_path(self):
        """The path to the database with the output results, None if not set."""
        try:
            return self._outdb_path
        except AttributeError:
            return None


class GbrvCompoundsFlow(Flow):
    """
    A Flow made of :class:`GbrvCompoundRelaxAndEosWork` works.
    Provides a finalize method that merges the work results and save
    them in "all_results.json".
    """
    def finalize(self):
        """
        This method is called when the flow is completed.
        Return 0 if success
        """
        # Read work results form JSON files and build list.
        all_results, dict_list = [], []
        for work in self:
            with open(work.outdir.path_in("gbrv_results.json"), "rt") as fh:
                results = json.load(fh)
                row = dict(
                    a0_abs_err=results["a0_abs_err"],
                    a0_rel_err=results["a0_rel_err"],
                )
                row.update({p.symbol: p.basename for p in work.pseudos})
                dict_list.append(row)
                all_results.append(results)

        # Build pandas DataFrame and print it
        frame = pd.DataFrame(dict_list)
        print(frame)
        #from tabulate import tabulate
        #print(tabulate(table, headers=["Pseudos", "a0_rel_err", "a0_abs_err"]))
        data = {"all_results": all_results, "frame": frame.to_json()}

        with open(self.outdir.path_in("all_results.json"), "wt") as fh:
            json.dump(data, fh) #, indent=-1, sort_keys=True)

        return super(GbrvCompoundsFlow, self).finalize()
