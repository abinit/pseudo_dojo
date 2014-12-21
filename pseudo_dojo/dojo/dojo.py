# coding: utf-8
from __future__ import division, print_function, unicode_literals

import os

from pymatgen.io.abinitio.tasks import TaskManager
from pymatgen.io.abinitio.flows import Flow
from pymatgen.io.abinitio.pseudos import Pseudo
from pseudo_dojo.dojo.works import PPConvergenceFactory, DeltaFactory, GbrvFactory

import logging
logger = logging.getLogger(__file__)

ALL_ACCURACIES = ("low", "normal", "high")

ALL_TRIALS = (
    "deltafactor",
    "gbrv_bcc",
    "gbrv_fcc",
)


class Dojo(object):
    """This object build the flows for the analysis/validation of pseudopotentials"""
    def __init__(self, workdir=None, manager=None, trials=ALL_TRIALS, accuracies=ALL_ACCURACIES):
        """
        Args:
            workdir: Working directory.
            manager: :class:`TaskManager` object that will handle the sumbmission of the jobs.
        """
        self.workdir = os.path.abspath(workdir) if workdir is not None else os.path.join(os.getcwd(), "DOJO")
        self.manager = TaskManager.from_user_config() if manager is None else manager
        self.trials, self.accuracies = trials, accuracies

        # List of pseudos analyzed by the Dojo and corresponding flows.
        self.pseudos, self.flows = [], []

    def add_pseudo(self, pseudo):
        """Add a pseudo to the Dojo."""
        pseudo = Pseudo.as_pseudo(pseudo)
        dojo_report = pseudo.dojo_report

        # Construct the flow 
        flow_workdir = os.path.join(self.workdir, pseudo.name)
        flow = Flow(workdir=flow_workdir, manager=self.manager)

        # Construct the flow according to the info found in the dojo report.
        if not pseudo.has_hints:
            # We need the hints in order to run the other tests
            factory = PPConvergenceFactory()
            ecut_work = factory.work_for_pseudo(pseudo, ecut_slice=slice(4, None, 1), nlaunch=4)
            flow.register_work(ecut_work)

        else:
            # Hints are available --> construct a flow for the different trials.
            dojo_trial = "deltafactor"
            if dojo_trial in self.trials:
                # Do we have this element in the deltafactor database?
                #if not df_database().has_symbol(pseudo.symbol):
                #    logger.warning("Cannot find %s in deltafactor database." % pseudo.symbol)

                delta_factory = DeltaFactory()
                kppa = 6750 # 6750 is the value used in the deltafactor code.
                kppa = 1

                for accuracy in self.accuracies:
                    if dojo_report.has_trial(dojo_trial, accuracy): continue
                    ecut, pawecutdg = self._ecut_pawecutdg(pseudo, accuracy)
                    work = delta_factory.work_for_pseudo(pseudo, accuracy=accuracy, kppa=kppa, ecut=ecut, pawecutdg=pawecutdg)

                    logger.info("Adding work for %s with accuracy %s" % (dojo_trial, accuracy))
                    work.set_dojo_accuracy(accuracy)
                    flow.register_work(work)

            # Test if GBRV tests are wanted.
            gbrv_structs = [s.split("_")[1] for s in self.trials if s.startswith("gbrv_")]

            if gbrv_structs:
                gbrv_factory = GbrvFactory()
                for struct_type in gbrv_structs:
                    dojo_trial = "gbrv_" + struct_type
                    for accuracy in self.accuracies:
                        if dojo_report.has_trial(dojo_trial, accuracy): continue
                        ecut, pawecutdg = self._ecut_pawecutdg(pseudo, accuracy)
                        work = gbrv_factory.relax_and_eos_work(pseudo, struct_type, ecut=ecut, pawecutdg=pawecutdg)

                        logger.info("Adding work for %s with accuracy %s" % (dojo_trial, accuracy))
                        work.set_dojo_accuracy(accuracy)
                        flow.register_work(work)

        flow.allocate()
        self.pseudos.append(pseudo)
        self.flows.append(flow)

    def _ecut_pawecutdg(self, pseudo, accuracy):
        hint = pseudo.hint_for_accuracy(accuracy=accuracy)
        return hint.ecut, hint.ecut * hint.aug_ratio

    def build(self):
        """Build the dojo."""
        for flow in self.flows:
            flow.build_and_pickle_dump()


class HintsAndGbrvDojo(Dojo):
    def __init__(self, workdir=None, manager=None):
        Dojo.__init__(self, workdir=workdir, manager=manager,
                      trials=("gbrv_bcc", "gbrv_fcc"), accuracies=["normal"])

