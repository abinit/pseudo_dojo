# coding: utf-8
"""Base class for Dojo Workflows."""
from __future__ import division, print_function, unicode_literals

import abc
import sys
import os
import numpy as np
from abipy import abilab

from monty.collections import AttrDict
from monty.pprint import pprint_table
from pymatgen.core.units import Ha_to_eV
from pymatgen.io.abinitio.eos import EOS
from pymatgen.io.abinitio.pseudos import Pseudo, PseudoTable
from pymatgen.io.abinitio.abiobjects import SpinMode, Smearing, KSampling, RelaxationMethod
from pymatgen.io.abinitio.works import Work, build_oneshot_phononwork, OneShotPhononWork
from abipy.core.structure import Structure
from pseudo_dojo.refdata.gbrv import gbrv_database
from pseudo_dojo.refdata.deltafactor import df_database, df_compute
from pseudo_dojo.dojo.gbrv_outdb import GbrvOutdb #, RocksaltOutdb GbrvRecord, 

import logging
logger = logging.getLogger(__name__)


class GbrvCompoundsFactory(object):
    """Factory class producing :class:`Work` objects for GBRV calculations."""
    def __init__(self):
        self._db = gbrv_database()

    def make_ref_structure(self, formula, struct_type, ref):
        """
        Return the structure used in the GBRV tests given the chemical formula, the structure type
        and the reference code.
        """
        # Get the entry in the database
        entry = self._db.get_entry(formula, struct_type)

        if entry is None: 
            logger.critical("Cannot find entry for %s, returning None!" % formula)
            return None

        # Build the structure and handle a possibly missing value.
        structure = entry.build_structure(ref=ref)

        if structure is None:
            logger.warning("No AE structure for %s\n Will use gbrv_uspp data." % formula)
            structure = entry.build_structure(ref="gbrv_uspp")
        
        if structure is None: 
            logger.critical("Cannot initialize structure for %s, returning None!" % formula)

        return structure

    def relax_and_eos_work(self, accuracy, pseudos, formula, struct_type, ecut=None, pawecutdg=None, ref="ae", **kwargs):
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
        pseudos = PseudoTable.as_table(pseudos)

        if pseudos.allpaw and pawecutdg is None:
            raise ValueError("pawecutdg must be specified for PAW calculations.")

        structure = self.make_ref_structure(formula, struct_type=struct_type, ref=ref)
 
        return GbrvCompoundRelaxAndEosWork(structure, formula, struct_type, pseudos, accuracy,
                                           ecut=ecut, pawecutdg=pawecutdg, **kwargs)


class GbrvCompoundRelaxAndEosWork(Work):

    def __init__(self, structure, formula, struct_type, pseudos, accuracy, ecut=None, pawecutdg=None, ngkpt=(8, 8, 8),
                 spin_mode="unpolarized", toldfe=1.e-9, smearing="fermi_dirac:0.001 Ha",
                 ecutsm=0.05, chksymbreak=0, workdir=None, manager=None, **kwargs):
        """
        Build a :class:`Work` for the computation of the relaxed lattice parameter.

        Args:   
            structure: :class:`Structure` object 
            structure_type: fcc, bcc 
            pseudos: Pseudopotentials
            ecut: Cutoff energy in Hartree
            ngkpt: MP divisions.
            spin_mode: Spin polarization mode.
            toldfe: Tolerance on the energy (Ha)
            smearing: Smearing technique.
            workdir: String specifing the working directory.
            manager: :class:`TaskManager` responsible for the submission of the tasks.
        """
        super(GbrvCompoundRelaxAndEosWork, self).__init__(workdir=workdir, manager=manager)

        self.pseudos = pseudos
        self.formula = formula
        self.struct_type = struct_type
        self.accuracy = accuracy
        self.set_name("_".join(["gbrv", struct_type, formula, accuracy]))

        # nband must be large enough to accomodate fractional occupancies.
        fband = kwargs.pop("fband", None)

        # FIXME
        #nband = gbrv_nband(self.pseudo)

        # TODO: toldfe for the EOS, tolvrs for the relaxation.
        # Set extra_abivars.
        self.extra_abivars = dict(
            ecut=ecut,
            pawecutdg=pawecutdg,
            toldfe=toldfe,
            #ecutsm=0.5,
            #nband=nband,
            paral_kgb=paral_kgb
        )
                                       
        self.extra_abivars.update(**kwargs)
        self.ecut = ecut
        self.smearing = Smearing.as_smearing(smearing)

        # Kpoint sampling: shiftk depends on struct_type
        #shiftk = {"fcc": [0, 0, 0], "bcc": [0.5, 0.5, 0.5]}.get(struct_type)
        shiftk = [0, 0, 0]
        #ngkpt = (4,4,4)

        self.ksampling = KSampling.monkhorst(ngkpt, chksymbreak=chksymbreak, shiftk=shiftk)
        self.spin_mode = SpinMode.as_spinmode(spin_mode)
        relax_algo = RelaxationMethod.atoms_and_cell()

        inp = abilab.AbinitInput(structure, pseudos)
        inp.add_abiobjects(self.ksampling, relax_algo, self.spin_mode, self.smearing)
        inp.set_vars(self.extra_abivars)

        # Register structure relaxation task.
        self.relax_task = self.register_relax_task(inp)

    def set_outdb(self, path):
        """
        This function set the outdb property (a database with the Gbrv results)
        Use this function when you want the work to write the results of the 
        calculation to the ouddb calculation.
        """
        self._outdb_path = path
        return self

    @property
    def outdb_path(self):
        """The database with the output results, None if not set."""
        try:
            return self._outdb_path
        except AttributeError:
            return None

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

            scf_input = abilab.AbinitInput(new_structure, self.pseudos)
            scf_input.add_abiobjects(self.ksampling, self.spin_mode, self.smearing)
            scf_input.set_vars(extra)

            # Register new task
            self.register_scf_task(scf_input)

        # Allocate new tasks and update the pickle database.
        self.flow.allocate()
        self.flow.build_and_pickle_dump()

    def compute_eos(self):
        self.history.info("Computing EOS")
        #results = self.get_results()

        # Read etotals and fit E(V) with a parabola to find the minimum
        etotals = self.read_etotals(unit="eV")[1:]
        assert len(etotals) == len(self.volumes)

        results = {}
        results.update(dict(
            etotals=list(etotals),
            volumes=list(self.volumes),
            num_sites=len(self.relaxed_structure),
        ))

        try:
            eos_fit = EOS.Quadratic().fit(self.volumes, etotals)
        except EOS.Error as exc:
            results.push_exceptions(exc)

        #return results

        # Function to compute cubic a0 from primitive v0 (depends on struct_type)
        vol2a = {"fcc": lambda vol: (4 * vol) ** (1/3.),
                 "rocksalt": lambda vol: (4 * vol) ** (1/3.),
                 "bcc": lambda vol: (2 * vol) ** (1/3.),
                 "ABO3": lambda vol: vol ** (1/3.),
                 }[self.struct_type]

        a0 = vol2a(eos_fit.v0)

        results.update(dict(
            v0=eos_fit.v0,
            b0=eos_fit.b0,
            b1=eos_fit.b1,
            a0=a0,
            ecut=self.ecut,
            #struct_type=self.struct_type
        ))

        # Update the database.
        # TODO, handle error!
        if self.outdb_path is not None:
            GbrvOutdb.update_record(self.outdb_path, self.formula, self.accuracy, self.pseudos, results)

        db = gbrv_database()
        entry = db.get_entry(self.formula, stype=self.struct_type)

        pawabs_err = a0 - entry.gbrv_paw
        pawrel_err = 100 * (a0 - entry.gbrv_paw) / entry.gbrv_paw

        # If AE results are missing we use GBRV_PAW as reference.
        if entry.ae is not None:
            abs_err = a0 - entry.ae
            rel_err = 100 * (a0 - entry.ae) / entry.ae
        else:
            abs_err = pawabs_err
            rel_err = pawrel_err

        print("for %s (%s) a0=%.2f Angstrom" % (self.formula, self.struct_type, a0))
        print("AE - THIS: abs_err = %f, rel_err = %f %%" % (abs_err, rel_err))
        print("GBRV-PAW - THIS: abs_err = %f, rel_err = %f %%" % (pawabs_err, pawrel_err))

        #d = {k: results[k] for k in ("a0", "etotals", "volumes")}
        #d["a0_abs_err"] = abs_err
        #d["a0_rel_err"] = rel_err
        #if results.exceptions:
        #    d["_exceptions"] = str(results.exceptions)

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
