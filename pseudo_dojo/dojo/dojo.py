from __future__ import division, print_function

import sys
import os
import numpy as np

from pymatgen.util.string_utils import pprint_table
from pymatgen.io.abinitio.eos import EOS
from pymatgen.io.abinitio.tasks import TaskManager
from pymatgen.io.abinitio.flows import AbinitFlow
from pymatgen.io.abinitio.pseudos import Pseudo, read_dojo_report
from pseudo_dojo.dojo.dojo_workflows import PPConvergenceFactory, DeltaFactory, GbrvFactory

import logging
logger = logging.getLogger(__file__)

ALL_ACCURACIES = ("low", "normal", "high")

ALL_TRIALS = (
    "deltafactor",
    "gbrv_bcc",
    "gbrv_fcc",
)


class DojoReport(dict):

    _TRIALS2KEY = {
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

    def has_exceptions(self):
        problems = {}

        for trial in self.ALL_TRIALS:
            for accuracy in ALL_ACCURACIES:
                excs = self[trial][accuracy].get("_exceptions", None)
                if excs is not None:
                    if trial not in problems:
                        problems[trial] = {}

                    problems[trial][accuracy] = excs

        return problems

    @property
    def has_hints(self):
        return "hints" in self

    def trials(self):
        return [k for k in self.keys() if  k != "hints"]

    def has_trial(self, dojo_trial, accuracy):
        """
        True if the dojo_report contains an entry for the
        given dojo_trial with the specified accuracy.
        If accuracy is None, we test if all accuracies are present
        """
        if dojo_trial not in self:
            return False

        if accuracy is not None:
            return accuracy in self[dojo_trial]
        else:
            return all(acc in self[dojo_trial] for acc in ALL_ACCURACIES)

    def to_table(self, **kwargs):
        """
        ===========  ===============  ===============   ===============
        Trial             low              normal            high 
        ===========  ===============  ===============   ===============
        deltafactor  value (rel_err)  value (rel_err)   value (rel_err)
        gbrv_fcc     ...              ...               ...
        ===========  ===============  ===============   ===============
        """
        # Build the header
        if kwargs.pop("with_hints", True):
            ecut_acc = {acc: self["hints"][acc]["ecut"] for acc in ALL_ACCURACIES}
            l = ["%s (%s Ha)" % (acc, ecut_acc[acc]) for acc in ALL_ACCURACIES]
        else:
            l = list(ALL_ACCURACIES)

        table = [["Trial"] + l]
        #row = ["%s (%s)" % (accuracy, ecut)]

        for trial in ALL_TRIALS:
            row = [trial]
            for accuracy in ALL_ACCURACIES:
                if not self.has_trial(trial, accuracy): 
                    row.append("N/A")
                else:
                    d = self[trial][accuracy]
                    #print(d.keys())
                    #s = "%s (%s %%)" % (value, rel_err)
                    value = d[self._TRIALS2KEY[trial]]
                    s = "%.1f" % value
                    row.append(s)

            table.append(row)

        return table

    def print_table(self, stream=sys.stdout):
        pprint_table(self.to_table(), out=stream)

    def plot_etotal_vs_ecut(self, **kwargs):
        d = self["hints"]
        ecuts, etotals, aug_ratios = np.array(d["ecuts"]), np.array(d["etotals"]), d["aug_ratio"]
        plot_etotals(ecuts, etotals, aug_ratios, **kwargs)

    def plot_gbrv_eos(self, struct_type, **kwargs):
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

        trial = "gbrv_" + struct_type
        for accuracy in ALL_ACCURACIES:
            if not self.has_trial(trial, accuracy): continue
            d = self[trial][accuracy]
            #num_sites, volumes, etotals = d["num_sites"], np.array(d["volumes"]), np.array(d["etotals"])
            volumes, etotals = np.array(d["volumes"]), np.array(d["etotals"])

            eos_fit = EOS.Quadratic().fit(volumes, etotals)
            eos_fit.plot(ax=ax, show=False) 

        plt.show()

    def plot_delta_factor_eos(self, **kwargs):
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

        trial = "deltafactor"
        for accuracy in ALL_ACCURACIES:
            if not self.has_trial(trial, accuracy): continue
            d = self[trial][accuracy]
            num_sites, volumes, etotals = d["num_sites"], np.array(d["volumes"]), np.array(d["etotals"])

            # Use same fit as the one employed for the deltafactor.
            eos_fit = EOS.DeltaFactor().fit(volumes/num_sites, etotals/num_sites)
            eos_fit.plot(ax=ax, show=False)

        plt.show()


class Dojo(object):
    """
    This object build the flows for the analysis/validation of pseudopotentials
    """
    def __init__(self, workdir=None, manager=None, trials=ALL_TRIALS, accuracies=ALL_ACCURACIES):
        """
        Args:
            workdir:
                Working directory.
            manager:
                `TaskManager` object that will handle the sumbmission of the jobs.
        """
        self.workdir = os.path.abspath(workdir) if workdir is not None else os.path.join(os.getcwd(), "DOJO")
        self.manager = TaskManager.from_user_config() if manager is None else manager
        self.trials, self.accuracies = trials, accuracies

        # List of pseudos analyzed by the Dojo and corresponding flows.
        self.pseudos, self.flows = [], []

    def add_pseudo(self, pseudo):
        """Add a pseudo to the Dojo."""
        pseudo = Pseudo.as_pseudo(pseudo)

        dojo_report = DojoReport.from_file(pseudo.filepath)

        # Construct the flow 
        flow_workdir = os.path.join(self.workdir, pseudo.name)
        flow = AbinitFlow(workdir=flow_workdir, manager=self.manager, pickle_protocol=0)

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
    for aratio, etot in zip(aug_ratios, etotals):
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
