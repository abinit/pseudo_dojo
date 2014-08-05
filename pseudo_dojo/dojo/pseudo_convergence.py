"""
Abinit Workflows
"""
from __future__ import division, print_function

import sys
import os
import shutil
import time
import abc
import collections
import numpy as np

from pymatgen.core.units import ArrayWithUnit, Ha_to_eV
from pymatgen.core.design_patterns import AttrDict
from pymatgen.serializers.json_coders import MSONable, json_pretty_dump
from pymatgen.util.num_utils import iterator_from_slice, chunks, monotonic
from pymatgen.util.string_utils import pprint_table, WildCard
from pymatgen.io.abinitio import wrappers
from pymatgen.io.abinitio.tasks import (Task, AbinitTask, Dependency, Node, ScfTask, NscfTask, BseTask, RelaxTask)
from pymatgen.io.abinitio.strategies import HtcStrategy, ScfStrategy #, RelaxStrategy
from pymatgen.io.abinitio.utils import Directory
from pymatgen.io.abinitio.netcdf import ETSF_Reader
from pymatgen.io.abinitio.abiobjects import Smearing, AbiStructure, KSampling, Electrons  #, RelaxationMethod
from pymatgen.io.abinitio.pseudos import Pseudo
from pymatgen.io.abinitio.abitimer import AbinitTimerParser
from pymatgen.io.abinitio.abiinspect import yaml_read_kpoints
from pymatgen.io.abinitio.workflows import Workflow


import logging
logger = logging.getLogger(__name__)


def check_conv(values, tol, min_numpts=1, mode="abs", vinf=None):
    """
    Given a list of values and a tolerance tol, returns the leftmost index for which

        abs(value[i] - vinf) < tol if mode == "abs"

    or

        abs(value[i] - vinf) / vinf < tol if mode == "rel"

    returns -1 if convergence is not achieved. By default, vinf = values[-1]

    Args:
        tol:
            Tolerance
        min_numpts:
            Minimum number of points that must be converged.
        mode:
            "abs" for absolute convergence, "rel" for relative convergence.
        vinf:
            Used to specify an alternative value instead of values[-1].
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

    ihigh   = check_conv(etotals, de_high, min_numpts=min_numpts)
    inormal = check_conv(etotals, de_normal)
    ilow    = check_conv(etotals, de_low)

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
        ecut_low    = ecuts[ilow]
        ecut_normal = ecuts[inormal]
        ecut_high   = ecuts[ihigh]

    aug_ratios = [1,]
    aug_ratio_low, aug_ratio_normal, aug_ratio_high = 3 * (1,)

    data = {
        "exit": ihigh != -1,
        "etotals": list(etotals),
        "ecuts": list(ecuts),
        "aug_ratios": aug_ratios,
        "low": {"ecut": ecut_low, "aug_ratio": aug_ratio_low},
        "normal": {"ecut": ecut_normal, "aug_ratio": aug_ratio_normal},
        "high": {"ecut": ecut_high, "aug_ratio": aug_ratio_high},
    }

    return AttrDict(**data)


def plot_etotals(ecuts, etotals, aug_ratios, **kwargs):
    """
    Uses Matplotlib to plot the energy curve as function of ecut

    Args:
        ecuts:
            List of cutoff energies
        etotals:
            Total energies in Hartree, see aug_ratios
        aug_ratios:
            List augmentation rations. [1,] for norm-conserving, [4, ...] for PAW
            The number of elements in aug_ration must equal the number of (sub)lists
            in etotals. Example:

                - NC: etotals = [3.4, 4,5 ...], aug_ratios = [1,]
                - PAW: etotals = [[3.4, ...], [3.6, ...]], aug_ratios = [4,6]

        =========     ==============================================================
        kwargs        description
        =========     ==============================================================
        show          True to show the figure
        savefig       'abc.png' or 'abc.eps'* to save the figure to a file.
        =========     ==============================================================

    Returns:
        `matplotlib` figure.
    """
    show = kwargs.pop("show", True)
    savefig = kwargs.pop("savefig", None)

    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    npts = len(ecuts)

    if len(aug_ratios) != 1 and len(aug_ratios) != len(etotals):
        raise ValueError("The number of sublists in etotal must equal the number of aug_ratios")

    if len(aug_ratios) == 1:
        etotals = [etotals,]

    lines, legends = [], []

    #emax = -np.inf
    for (aratio, etot) in zip(aug_ratios, etotals):
        emev = np.array(etot) * Ha_to_eV * 1000
        emev_inf = npts * [emev[-1]]
        yy = emev - emev_inf
        #print("emax", emax)
        #print("yy", yy)
        #emax = np.max(emax, np.max(yy))

        line, = ax.plot(ecuts, yy, "-->", linewidth=3.0, markersize=10)

        lines.append(line)
        legends.append("aug_ratio = %s" % aratio)

    ax.legend(lines, legends, 'upper right', shadow=True)

    # Set xticks and labels.
    ax.grid(True)
    ax.set_title("$\Delta$ Etotal Vs Ecut")
    ax.set_xlabel("Ecut [Ha]")
    ax.set_ylabel("$\Delta$ Etotal [meV]")
    ax.set_xticks(ecuts)

    #ax.yaxis.set_view_interval(-10, emax + 0.01 * abs(emax))
    ax.yaxis.set_view_interval(-10, 20)

    if show:
        plt.show()

    if savefig is not None:
        fig.savefig(savefig)

    return fig


class PseudoConvergence(Workflow):

    def __init__(self, pseudo, ecuts, atols_mev,
                 toldfe=1.e-8, spin_mode="polarized", acell=(8, 9, 10), 
                 smearing="fermi_dirac:0.1 eV", workdir=None, manager=None):

        super(PseudoConvergence, self).__init__(workdir, manager)

        self.pseudo = Pseudo.as_pseudo(pseudo)
        self.ecuts = list(ecuts)
        self.atols_mev = atols_mev
        self.toldfe = toldfe
        self.spin_mode = spin_mode
        self.acell = acell
        self.smearing = smearing

        for ecut in self.ecuts:
            self.add_task_with_ecut(ecut)

    def add_task_with_ecut(self, ecut):
        """Register a new task with cutoff energy ecut."""
        # One atom in a box of lenghts acell.
        boxed_atom = AbiStructure.boxed_atom(self.pseudo, acell=self.acell)

        # Gamma-only sampling.
        gamma_only = KSampling.gamma_only()

        # Don't write WFK files.
        extra_abivars = {
            "ecut" : ecut,
            "prtwf": 0,
            "toldfe": self.toldfe,
        }

        strategy = ScfStrategy(boxed_atom, self.pseudo, gamma_only,
                               spin_mode=self.spin_mode, smearing=self.smearing,
                               **extra_abivars)

        self.register(strategy)

    def get_results(self):
        # Get the results of the tasks.
        wf_results = super(PseudoConvergence, self).get_results()

        etotals = self.read_etotals()
        data = compute_hints(self.ecuts, etotals, self.atols_mev)
        #print("ecuts", self.ecuts)
        #print("etotals", etotals)

        #plot_etotals(data["ecuts"], data["etotals"], data["aug_ratios"],
        #            show=False, savefig=self.path_in_workdir("etotals.pdf"))

        wf_results.update(data)

        if not monotonic(etotals, mode="<", atol=1.0e-5):
            logger.warning("E(ecut) is not decreasing")
            wf_results.push_exceptions("E(ecut) is not decreasing:\n" + str(etotals))

        return wf_results


class PseudoIterativeConvergence(Workflow):

    def __init__(self, pseudo, ecut_slice, atols_mev,
                 toldfe=1.e-8, spin_mode="polarized", acell=(8, 9, 10), 
                 smearing="fermi_dirac:0.1 eV", max_niter=50, workdir=None, manager=None):
        """
        Args:
            pseudo:
                string or Pseudo instance
            ecut_slice:
                List of cutoff energies or slice object (mainly used for infinite iterations).
            atols_mev:
                List of absolute tolerances in meV (3 entries corresponding to accuracy ["low", "normal", "high"]
            spin_mode:
                Defined how the electronic spin will be treated.
            acell:
                Lengths of the periodic box in Bohr.
            smearing:
                Smearing instance or string in the form "mode:tsmear". Default: FemiDirac with T=0.1 eV
            workdir:
                Working directory.
            manager:
                `TaskManager` object.
        """
        super(PseudoIterativeConvergence, self).__init__(workdir, manager)

        self.pseudo = Pseudo.as_pseudo(pseudo)
        self.atols_mev = atols_mev
        self.toldfe = toldfe
        self.spin_mode = spin_mode
        self.acell = acell
        self.smearing = smearing
        self.ecuts = []

        assert isinstance(ecut_slice, slice)
        self.ecut_slice = ecut_slice

        if self.pseudo.ispaw:
            raise NotImplementedError("PAW convergence tests are not supported yet")

        self.add_task_with_ecut(ecut_slice.start)

    def add_task_with_ecut(self, ecut):
        """Register a new task with cutoff energy ecut."""
        # One atom in a box of lenghts acell.
        boxed_atom = AbiStructure.boxed_atom(self.pseudo, acell=self.acell)

        # Gamma-only sampling.
        gamma_only = KSampling.gamma_only()

        # Don't write WFK files.
        extra_abivars = {
            "ecut" : ecut,
            "prtwf": 0,
            "toldfe": self.toldfe,
        }

        strategy = ScfStrategy(boxed_atom, self.pseudo, gamma_only,
                               spin_mode=self.spin_mode, smearing=self.smearing,
                               **extra_abivars)

        self.ecuts.append(ecut)
        self.register(strategy)

    #def get_results(self):
    #    """Return the results of the tasks."""
    #    wf_results = super(PseudoIterativeConvergence, self).get_results()

    #    data = self.check_etotal_convergence()

    #    ecuts, etotals, aug_ratios = data["ecuts"],  data["etotals"], data["aug_ratios"]

    #    #plot_etotals(ecuts, etotals, aug_ratios, show=False, savefig=self.path_in_workdir("etotals.pdf"))

    #    wf_results.update(data)

    #    if not monotonic(data["etotals"], mode="<", atol=1.0e-5):
    #        logger.warning("E(ecut) is not decreasing")
    #        wf_results.push_exceptions("E(ecut) is not decreasing\n" + str(etotals))

    #    return wf_results

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
            etotals=results["etotals"],
        ))
        
        if results.exceptions:
            d["_exceptions"] = str(results.exceptions)

        return {self.dojo_key: d}

    def on_all_ok(self):
        """
        This method is called when self reaches S_OK.
        """
        etotals = self.read_etotals()
        data = compute_hints(self.ecuts, etotals, self.atols_mev)

        if data.exit:
            print("converged")

            #self.write_dojo_report(d)

            #elow, ehigh = data.low["ecut"], data.high["ecut"]
            #for ecut in range(elow, ehigh+1, 1):
            #    self.add_task_with_ecut(ecut)

            #self._finalized = False
            #self.flow.allocate()
            #self.flow.build_and_pickle_dump()
        else:
            print("Building new tasks")
            self.add_task_with_ecut(self.ecuts[-1] + self.ecut_slice.step)
            self._finalized = False
            self.flow.allocate()
            self.flow.build_and_pickle_dump()

        return super(PseudoIterativeConvergence, self).on_all_ok()


class PPConvergenceFactory(object):
    """
    Factory object that constructs workflows for analyzing the converge of
    pseudopotentials.
    """
    def work_for_pseudo(self, pseudo, ecut_slice, 
                        toldfe=1.e-8, atols_mev=(10, 1, 0.1), spin_mode="polarized",
                        acell=(8, 9, 10), smearing="fermi_dirac:0.1 eV", workdir=None, manager=None):
        """
        Return a `Workflow` object given the pseudopotential pseudo.

        Args:
            pseudo:
                Pseudo object.
            ecut_slice:
                cutoff energies in Ha units (accepts lists or slice objects)
            toldfe:
                Tolerance on the total energy (Ha).
            atols_mev:
                Tolerances in meV for accuracy in ["low", "normal", "high"]
            spin_mode:
                Spin polarization.
            acell:
                Length of the real space lattice (Bohr units)
            smearing:
                Smearing technique.
            workdir:
                Working directory.
            manager:
                `TaskManager` object.
        """
        if ecut_slice.stop is None:
            return PseudoIterativeConvergence(
                    pseudo, ecut_slice, atols_mev, 
                    toldfe=toldfe, spin_mode=spin_mode, 
                    acell=acell, smearing=smearing, workdir=workdir, manager=manager)

        else:
            ecuts = np.arange(ecut_slice.start, ecut_slice.stop, ecut_slice.step)
            return PseudoConvergence(
                        pseudo, ecuts, atols_mev, 
                        toldfe=toldfe, spin_mode=spin_mode, 
                        acell=acell, smearing=smearing, workdir=workdir, manager=manager)


#class HintsMaster(object): 
#    """
#    Level 0 master that analyzes the convergence of the total energy versus
#    the plane-wave cutoff energy.
#    """
#    dojo_level = 0
#    dojo_key = "hints"
#
#    # Absolute tolerance for low, normal, and high accuracy.
#    _ATOLS_MEV = (10, 1, 0.1)
#
#    def challenge(self, workdir):
#        workdir = os.path.join(workdir, "LEVEL_" + str(self.dojo_level))
#        flow = abilab.AbinitFlow(workdir=workdir, manager=self.manager)
#
#        from pseudo_dojo.dojo.pseudo_convergence import PPConvergenceFactory
#        factory = PPConvergenceFactory()
#
#        pseudo = self.pseudo
#        toldfe = 1.e-10
#        estep = 10 #kwargs.get("estep", 10)
#        eslice = slice(5, None, estep)
#
#        # TODO Rewrite this
#        #w = factory.work_for_pseudo(pseudo, eslice, toldfe=toldfe, atols_mev=self._ATOLS_MEV)
#
#        #estart = max(wres["low"]["ecut"] - estep, 5)
#        #if estart <= 10:
#        #    estart = 1 # To be sure we don't overestimate ecut_low
#
#        #estop, estep = wres["high"]["ecut"] + estep, 1
#        estart=10; estop=25; estep=5
#
#        erange = list(np.arange(estart, estop, estep))
#
#        self.work = work = factory.work_for_pseudo(pseudo, erange, toldfe=toldfe, atols_mev=self._ATOLS_MEV)
#
#        flow.register_work(work)
#        flow.allocate()
#        flow.build_and_pickle_dump()
#
#        print("Finding optimal values for ecut in the range [%.1f, %.1f, %1.f,] Hartree, "
#              % (estart, estop, estep))


def build_flow():
    from abipy import abilab
    factory = PPConvergenceFactory()
    pseudo = sys.argv[1]

    flow = abilab.AbinitFlow("hello_pseudo", abilab.TaskManager.from_user_config(), pickle_protocol=0)

    atols_mev = (10, 1, 0.1)
    ecut_slice = slice(5, None, 1)

    work = PseudoIterativeConvergence(pseudo, ecut_slice, atols_mev)
    #work = factory.work_for_pseudo(pseudo, [2, 4, 6])

    flow.register_work(work)

    return flow.allocate()
