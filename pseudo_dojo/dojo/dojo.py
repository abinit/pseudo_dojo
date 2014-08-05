from __future__ import division, print_function

import sys
import os
import numpy as np

from abipy import abilab
from pymatgen.util.string_utils import pprint_table
from pymatgen.io.abinitio.tasks import TaskManager
from pymatgen.io.abinitio.pseudos import Pseudo, read_dojo_report
from pseudo_dojo.dojo.dojo_workflows import PPConvergenceFactory, DeltaFactory, GbrvFactory

import logging
logger = logging.getLogger(__file__)

_ALL_ACCURACIES = ["low", "normal", "high"]


#class DojoReportSection(dict)
#    #dojo_trial
#
#    def __init__(self, dojo_accuracy, **kwargs)
#        super(ReportEntry, self).__init__(**kwargs)
#        self.dojo_accuracy = dojo_accuracy
#
#    def plot(self, ax=None, **kwargs)
#
#class DeltaFactorSection(DojoReportSection)
#    dojo_trial = "delta_factor"


class DojoReport(dict):
    _TESTS = [
        "deltafactor",
        "gbrv_bcc",
        "gbrv_fcc",

    ]

    _TESTS2KEY = {
        "deltafactor": "dfact_meV",
        "gbrv_bcc": "a0",
        "gbrv_fcc": "a0",
    }

    def __init__(self, filepath): 
        super(dict, self).__init__()
        d = read_dojo_report(filepath)
        self.update(**d)

    @classmethod
    def from_file(cls, filepath):
        return cls(filepath)

    #def groupby(self):
    #def has_trial(self, dojo_trial, accuracy=None)

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
                print(d.keys())
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


class Dojo(object):
    """
    This object drives the execution of the tests for the pseudopotential.

    A Dojo has a set of masters, each master is associated to a particular trial
    and the master is responsible for the validation/rating of the results of the tests.
    """

    ACCURACIES = [
        "low", 
        #"normal", 
        #"high"
    ]

    TRIALS = [
        #"deltafactor",
        "GbrvRelax",
    ]

    def __init__(self, workdir=None, manager=None):
        """
        Args:
            pseudo:
            workdir:
                Working directory.
            manager:
                `TaskManager` object that will handle the sumbmission of the jobs.
        """
        self.workdir = os.path.abspath(workdir) if workdir is not None else os.path.join(os.getcwd(), "DOJO")
        self.manager = TaskManager.from_user_config() if manager is None else manager

        # List of pseudos analyzed by the Dojo and corresponding flows.
        self.pseudos, self.flows = [], []

    #def __str__(self):

    def add_pseudo(self, pseudo):
        """Add a pseudo to the Dojo."""
        pseudo = Pseudo.as_pseudo(pseudo)
        #dojo_report = DojoReport(pseudo.filepath)
        #print(dojo_report)

        # Construct the flow 
        flow_workdir = os.path.join(self.workdir, pseudo.name)
        flow = abilab.AbinitFlow(workdir=flow_workdir, manager=self.manager, pickle_protocol=0)

        # Inspect the dojo_report and build the flow
        if not pseudo.has_hints:
            factory = PPConvergenceFactory()
            ecut_work = factory.work_for_pseudo(pseudo, ecut_slice=slice(4, None, 1), nlaunch=4)
            flow.register_work(ecut_work)

        else:
            hint = pseudo.hint_for_accuracy(accuracy="normal")
            ecut = hint.ecut
            pawecutdg = ecut * hint.aug_ratio

            dojo_trial = "deltafactor"
            if dojo_trial in self.TRIALS:
                # Do we have this element in the deltafactor database?
                #if not df_database().has_symbol(pseudo.symbol):
                #    logger.warning("Cannot find %s in deltafactor database." % pseudo.symbol)

                delta_factory = DeltaFactory()
                kppa = 6750 # 6750 is the value used in the deltafactor code.
                kppa = 1

                for accuracy in self.ACCURACIES:
                    #if dojo_report.has_trial(dojo_trial, accuracy=accuracy): continue
                    work = delta_factory.work_for_pseudo(pseudo, accuracy=accuracy, kppa=kppa, ecut=ecut, pawecutdg=pawecutdg)

                    work.set_dojo_accuracy(accuracy)
                    flow.register_work(work)

            if "GbrvRelax" in self.TRIALS:
                gbrv_factory = GbrvFactory()
                for struct_type in ["fcc", "bcc"]:
                    for accuracy in self.ACCURACIES:
                        #if dojo_report.has_trial(dojo_trial, accuracy=accuracy): continue
                        work = gbrv_factory.relax_and_eos_work(pseudo, struct_type, ecut=ecut, pawecutdg=pawecutdg)

                        work.set_dojo_accuracy(accuracy)
                        flow.register_work(work)

        flow.allocate()
        self.pseudos.append(pseudo)
        self.flows.append(flow)

    def build(self):
        """Build the dojo."""
        for flow in self.flows:
            flow.build_and_pickle_dump()
