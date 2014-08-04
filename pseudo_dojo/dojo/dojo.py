from __future__ import division, print_function

import sys
import os
import abc
import shutil
import numpy as np
from abipy import abilab

from pymatgen.util.string_utils import pprint_table
from pymatgen.io.abinitio.tasks import TaskManager
from pymatgen.io.abinitio.pseudos import Pseudo, read_dojo_report
from pymatgen.io.abinitio.calculations import PPConvergenceFactory
from pymatgen.io.abinitio.launcher import PyFlowScheduler
from pseudo_dojo.dojo.deltaworks import DeltaFactory
from pseudo_dojo.refdata.deltafactor import df_database, df_compute

import logging
logger = logging.getLogger(__file__)

_ALL_ACCURACIES = ["low", "normal", "high"]


class DojoReport(dict):
    _TESTS = ["delta_factor",]

    _TESTS2KEY = dict(
        delta_factor="dfact_meV",
    )

    def __init__(self, filepath): 
        super(dict, self).__init__()
        d = read_dojo_report(filepath)
        self.update(**d)

    def print_table(self, stream=sys.stdout):
        pprint_table(self.to_table(), out=stream)

    def to_table(self):
        """
        ========  ================ =========
        Accuracy  DeltaFactor      Gbrv_fcc
        ========  ================ =========
        low       value (rel_err)
        normal
        high 
        ========  ================ =========
        """
        #symbol = self["symbol"]
        table = [["Accuracy"] + self._TESTS]

        for accuracy in _ALL_ACCURACIES:
            #ecut = self["hints"][accuracy]["ecut"]
            #row = ["%s (%s)" % (accuracy, ecut)]
            row = [accuracy]
            for test in self._TESTS:
                d = self[test][accuracy]
                #value = d["dojo_value"]
                #rel_err = d["dojo_rel_err"]
                #s = "%s (%s %%)" % (value, rel_err)
                value = d[self._TESTS2KEY[test]]
                s = "%.1f" % value
                row.append(s)

            table.append(row)

        return table

    def find_exceptions(self):
        problems = {}

        for test in self._TESTS:
            for accuracy in _ALL_ACCURACIES:
                excs = self[test][accuracy].get("_exceptions", None)
                if excs is not None:
                    if test not in problems:
                        problems[test] = {}
                    problems[tests][accuracy] = excs

        return problems

    #def plot_etotal_vs_ecut(self, **kwargs):

    def plot_delta_factor_eos(self, **kwargs):
        # Use same fit as the one employed for the deltafactor.
        from pymatgen.io.abinitio.eos import EOS

        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

        for accuracy in _ALL_ACCURACIES:
            d = self["delta_factor"][accuracy]
            num_sites = d["num_sites"]
            volumes = np.array(d["volumes"])
            etotals = np.array(d["etotals"])
            eos_fit = EOS.DeltaFactor().fit(volumes/num_sites, etotals/num_sites)
            eos_fit.plot(ax=ax, show=False) #, savefig=self.outdir.path_in("eos.pdf"))

        plt.show()

    #def plot_gbrv_eos(self, **kwargs):


class DojoError(Exception):
    """Base Error class for DOJO calculations."""


class Dojo(object):
    """
    This object drives the execution of the tests for the pseudopotential.

    A Dojo has a set of masters, each master is associated to a particular trial
    and the master is responsible for the validation/rating of the results of the tests.
    """
    Error = DojoError

    def __init__(self, pseudo, workdir=None, manager=None, max_level=None):
        """
        Args:
            pseudo:
            workdir:
                Working directory.
            manager:
                `TaskManager` object that will handle the sumbmission of the jobs.
            max_level:
                Max test level to perform.
        """
        self.pseudo = Pseudo.as_pseudo(pseudo)
        #basedir = os.path.abspath(workdir)
        #self.workdir = os.path.abspath(workdir)
        self.workdir = os.path.join("DOJO_", self.pseudo.name)
        self.manager = TaskManager.from_user_config() if manager is None else manager

        # List of master classes that will be instantiated afterwards.
        # They are ordered according to the master level.
        classes = [m for m in DojoMaster.__subclasses__()]
        classes.sort(key=lambda cls: cls.dojo_level)

        self.master_classes = classes
        if max_level is not None:
            self.master_classes = classes[:max_level+1]

        #self.flows, self.masters = [], []
        #for master in all_masters:
        #    if master.accept_pseudo(self.pseudo):
        #        flow = master.build_flow(self.workdir)
        #        self.flows.append(flow)
        #        self.masters.append(master)

        #self.scheduler = PyFlowScheduler.from_user_config()
        #for flow in self.flows:
        #    self.scheduler.add_flow(flow)
        #self.scheduler.start()

    def start_training(self):
        """
        This method represents the main entry point for client code.
        The Dojo receives a pseudo-like object and delegate the execution
        of the tests to the dojo_masters
        """
        # Build master instances.
        all_masters = [cls(manager=self.manager) for cls in self.master_classes]

        isok = False
        for master in all_masters:
            if master.accept_pseudo(self.pseudo):
                isok = master.start_training(self.workdir)
                if not isok:
                    print("master: %s returned isok %s.\n Skipping next trials!" % (master.name, isok))
                    break
        return isok

    def __str__(self):
        return repr_dojo_levels()


class DojoMaster(object):
    """"
    Abstract base class for the dojo masters.
    Subclasses must define the class attribute level.
    """
    __metaclass__ = abc.ABCMeta

    Error = DojoError

    def __init__(self, manager):
        """
        Args:
            manager:
                `TaskManager` object 
        """
        self.manager = manager

    @property
    def name(self):
        """Name of the subclass."""
        return self.__class__.__name__

    @staticmethod
    def subclass_from_dojo_level(dojo_level):
        """Returns a subclass of `DojoMaster` given the dojo_level."""
        classes = []
        for cls in DojoMaster.__subclasses__():
            if cls.dojo_level == dojo_level:
                classes.append(cls)

        if len(classes) != 1:
            raise DojoError("Found %d masters with dojo_level %d" % (len(classes), dojo_level))

        return classes[0]

    def inspect_pseudo(self, pseudo):
        """Returns the maximum level of the DOJO trials passed by the pseudo."""
        pseudo.read_dojo_report()

        if not pseudo.has_dojo_report:
            max_level = None
        else:
            levels = [dojo_key2level(key) for key in pseudo.dojo_report]
            max_level = max(levels)

        return max_level

    def accept_pseudo(self, pseudo):
        """
        Returns True if the master can train the pseudo.
        This method is called before testing the pseudo.

        A master can train the pseudo if his level == pseudo.dojo_level + 1
        """
        ready = False
        pseudo_dojo_level = self.inspect_pseudo(pseudo)
        
        if pseudo_dojo_level is None:
            # Hints are missing
            ready = (self.dojo_level == 0)
            pseudo.write_dojo_report(report={})
        else:
            if pseudo_dojo_level == self.dojo_level:
                ready = False
                # pseudo has already a test associated to this level.
                # check if it has the same accuracy.
                #accuracy = kwargs.get("accuracy", "normal")
                #if accuracy not in pseudo.dojo_report[self.dojo_key]:
                #    ready = True
                #else:
                #    print("%s: %s has already an entry for accuracy %s" % (self.name, pseudo.name, accuracy))
                #    ready = False

            else:
                # Pseudo level must be one less than the level of the master.
                ready = (pseudo_dojo_level == self.dojo_level - 1)

        if not ready:
            print("%s: Sorry, %s-san, I cannot train you" % (self.name, pseudo.name))
        else:
            print("%s: Welcome %s-san, I'm your level-%d trainer" % (self.name, pseudo.name, self.dojo_level))
            self.pseudo = pseudo

        return ready

    def write_dojo_report(self, report, overwrite_data=False, ignore_errors=False):
        """
        Write/update the DOJO_REPORT section of the pseudopotential.
        """
        # Read old_report from pseudo.
        old_report = self.pseudo.read_dojo_report()

        dojo_key = self.dojo_key
        if dojo_key not in old_report:
            # Create new entry
            old_report[dojo_key] = {}
        else:
            pass
            # Check that we are not going to overwrite data.
            #if self.dojo_accuracy in old_report[dojo_key] and not overwrite_data:
            #    raise self.Error("%s already exists in the old pseudo. Cannot overwrite data" % dojo_key)

        # Update old report card with the new one.
        old_report[dojo_key].update(report[dojo_key])

        # Write new report
        self.pseudo.write_dojo_report(old_report)

    def start_training(self, workdir):
        """Start the tests in the working directory workdir."""
        self.challenge(workdir)
        report = self.make_report()
        self.write_dojo_report(report)

        isok = True
        if "_exceptions" in report:
            isok = False
            logger.warning("Found exceptions:\n %s" % str(report["_exceptions"]))

        return isok

    @abc.abstractmethod
    def challenge(self, workdir):
        """Abstract method to run the calculation."""

    @abc.abstractmethod
    def make_report(self):
        """
        Abstract method.
        Returns: 
            report:
                Dictionary with the results of the trial.
        """


class HintsMaster(DojoMaster):
    """
    Level 0 master that analyzes the convergence of the total energy versus
    the plane-wave cutoff energy.
    """
    dojo_level = 0
    dojo_key = "hints"

    # Absolute tolerance for low, normal, and high accuracy.
    _ATOLS_MEV = (10, 1, 0.1)

    def challenge(self, workdir):
        workdir = os.path.join(workdir, "LEVEL_" + str(self.dojo_level))
        flow = abilab.AbinitFlow(workdir=workdir, manager=self.manager)

        factory = PPConvergenceFactory()

        pseudo = self.pseudo
        toldfe = 1.e-10
        estep = 10 #kwargs.get("estep", 10)
        eslice = slice(5, None, estep)

        # TODO Rewrite this
        #w = factory.work_for_pseudo(pseudo, eslice, toldfe=toldfe, atols_mev=self._ATOLS_MEV)
        #flow.register_work(w)
        #flow.allocate()
        #flow.build_and_pickle_dump()
        #w = flow[0]
        #scheduler = PyFlowScheduler.from_user_config()
        #scheduler.add_flow(flow)

        #if os.path.exists(w.workdir): shutil.rmtree(w.workdir)
        #print("Converging %s in iterative mode with ecut_slice %s" % (pseudo.name, eslice, self.max_ncpus))
        #scheduler.start()

        #w.start()
        #w.wait()
        #wres = w.get_results()
        #w.move("ITERATIVE")

        #estart = max(wres["low"]["ecut"] - estep, 5)
        #if estart <= 10:
        #    estart = 1 # To be sure we don't overestimate ecut_low

        #estop, estep = wres["high"]["ecut"] + estep, 1
        estart=10; estop=25; estep=5

        erange = list(np.arange(estart, estop, estep))

        self.work = work = factory.work_for_pseudo(pseudo, erange, toldfe=toldfe, atols_mev=self._ATOLS_MEV)

        flow.register_work(work)
        flow.allocate()
        flow.build_and_pickle_dump()

        print("Finding optimal values for ecut in the range [%.1f, %.1f, %1.f,] Hartree, "
              % (estart, estop, estep))

        scheduler = PyFlowScheduler.from_user_config()
        scheduler.add_flow(flow)
        scheduler.start()

    def make_report(self):
        """
        "hints": {
            "high": {"aug_ratio": 1, "ecut": 45},
            "low": {...},
            "normal": {...}
        """
        results = self.work.get_results()
        d = {key: results[key] for key in _ALL_ACCURACIES}

        d.update(dict(
            ecuts=results["ecuts"],
            etotals=results["etotals"],
        ))
        
        if results.exceptions:
            d["_exceptions"] = str(results.exceptions)

        return {self.dojo_key: d}


class DeltaFactorMaster(DojoMaster):
    """
    This master performs the computation of the delta factor.
    """
    dojo_level = 1
    dojo_key = "delta_factor"

    def accept_pseudo(self, pseudo):
        """Returns True if the master can train the pseudo."""
        ready = super(DeltaFactorMaster, self).accept_pseudo(pseudo)

        # Do we have this element in the deltafactor database?
        return (ready and df_database().has_symbol(self.pseudo.symbol))

    def challenge(self, workdir):
        # Calculations will be executed in this directory.
        workdir = os.path.join(workdir, self.dojo_key)
        self.flow = flow = abilab.AbinitFlow(workdir=workdir, manager=self.manager, pickle_protocol=0)

        # 6750 is the value used in the deltafactor code.
        kppa = 6750
        kppa = 1

        factory = DeltaFactory()
        for accuracy in _ALL_ACCURACIES:
            work = factory.work_for_pseudo(self.pseudo, accuracy=accuracy, kppa=kppa, ecut=None, pawecutdg=None)
            work.dojo_accuracy = accuracy
            flow.register_work(work)

        flow.allocate()
        flow.build_and_pickle_dump()
        #scheduler = PyFlowScheduler.from_user_config()
        #scheduler.add_flow(flow)
        #scheduler.start()

    def make_report(self):
        # Get reference results (Wien2K).
        wien2k = df_database().get_entry(self.pseudo.symbol)

        data = {}
        for work in self.flow:
            accuracy = work.dojo_accuracy
            results = work.get_results()
            #results.json_dump(work.path_in_workdir("dojo_results.json"))

            # Get our results and compute deltafactor estimator.
            #v0, b0_GPa, b1 = results["v0"], results["b0_GPa"], results["b1"]
            #dfact = df_compute(wien2k.v0, wien2k.b0_GPa, wien2k.b1, v0, b0_GPa, b1, b0_GPa=True)
            #print("Deltafactor = %.3f meV" % dfact)

            d = {k: results[k] for k in ("dfact_meV", "v0", "b0", "b0_GPa", "b1", "etotals", "volumes", "num_sites")}

            if results.exceptions:
                d["_exceptions"] = str(results.exceptions)

            data.update({accuracy: d})

        return {self.dojo_key: data}


_key2level = {}
for cls in DojoMaster.__subclasses__():
    _key2level[cls.dojo_key] = cls.dojo_level


def dojo_key2level(key):
    """Return the trial level from the name found in the pseudo."""
    return _key2level[key]


def repr_dojo_levels():
    """String representation of the different levels of the Dojo."""
    level2key = {v: k for k, v in _key2level.items()}

    lines = ["Dojo level --> Challenge"]
    for k in sorted(level2key):
        lines.append("level %d --> %s" % (k, level2key[k]))

    return "\n".join(lines)
