"""
"""
from __future__ import unicode_literals, division, print_function

import sys
import json
import logging
import copy
import numpy as np

from collections import OrderedDict, defaultdict, Iterable
from tabulate import tabulate
from monty.string import list_strings, is_string
from monty.bisect import find_le
from pymatgen.analysis.eos import EOS
from pymatgen.core.periodic_table import Element
from pymatgen.util.plotting_utils import add_fig_kwargs, get_ax_fig_plt

logger = logging.getLogger(__name__)


class DojoReportError(Exception):
    """Exception raised by DoJoReport."""


class DojoReport(dict):
    """ 
    Dict-like object with the validation results.
    """

    _TRIALS2KEY = {
        "deltafactor": "dfact_meV",
        "gbrv_bcc": "a0_rel_err",
        "gbrv_fcc": "a0_rel_err",
        "phwoa": "all",
        "phonon": "all",
        "ebands": "all"
    }

    ALL_ACCURACIES = ("low", "normal", "high")

    # List of dojo_trials
    # Remember to update the list if you add a new test to the DOJO_REPORT
    ALL_TRIALS = (
        "deltafactor",
        "gbrv_bcc",
        "gbrv_fcc",
        "phonon",
        "phwoa",
        "ebands"
    )

    # Tolerances on the deltafactor prime (in eV) used for the hints.
    ATOLS = (0.5, 0.1, 0.02)
    #for noble gasses:
    #ATOLS = (1.0, 0.2, 0.04)

    LAST_VERSION = "1.0"

    Error = DojoReportError

    @classmethod
    def from_file(cls, filepath):
        """Read the DojoReport from file."""
        with open(filepath, "rt") as fh:
            return cls(**json.load(fh))

    @classmethod
    def new_from_pseudo(cls, pseudo):
        """Build a DojoReport from Pseudo."""
        new = cls()
        new["symbol"] = pseudo.symbol
        new["md5"]  = pseudo.compute_md5()
        new["version"] = cls.LAST_VERSION
        return new

    @classmethod
    def from_hints(cls, ppgen_ecut, symbol):
        """
        Initialize an empty DojoReport from the initial guesses for  
        the cutoff energy in Hartree

        Args:
            ppgen_ecut: tuple(3) cutoff energies for the 3 accuracy levels.
            symbol: Chemical symbol.
        """
        dense_right = np.arange(ppgen_ecut, ppgen_ecut + 6*2, step=2)
        dense_left = np.arange(max(ppgen_ecut-6, 2), ppgen_ecut, step=2)
        coarse_high = np.arange(ppgen_ecut + 15, ppgen_ecut + 35, step=5)

        ecuts = list(dense_left) + list(dense_right) + list(coarse_high)
        return cls(ecuts=ecuts, symbol=symbol) 

    def __init__(self, *args, **kwargs):
        super(DojoReport, self).__init__(*args, **kwargs)

        try:
            for trial in self.ALL_TRIALS:
                # Convert ecut to float and build an OrderedDict (results are indexed by ecut in ascending order)
                try:
                    d = self[trial]
                except KeyError:
                    continue
                ecuts_keys = sorted([(float(k), k) for k in d], key=lambda t: t[0])
                self[trial] = OrderedDict([(t[0], d[t[1]]) for t in ecuts_keys])

        except ValueError:
            raise self.Error('Error while initializing the dojo report')

        if "version" not in self:
            self["version"] = self.LAST_VERSION

    def __str__(self):
        return(json.dumps(self, indent=-1))

    def deepcopy(self):
        """Deepcopy of the object."""
        return copy.deepcopy(self)

    @property
    def symbol(self):
        """Chemical symbol."""
        return self["symbol"]

    @property
    def element(self):
        """Element object."""
        return Element(self.symbol)

    @property
    def ecuts(self):
        """Numpy array with the list of ecuts that should be present in the dojo_trial sub-dicts"""
        return self["ecuts"]

    @property
    def trials(self):
        """Set of strings with the trials present in the report."""
        return set(list(self.keys())).intersection(self.ALL_TRIALS)

    def has_trial(self, dojo_trial, ecut=None):
        """
        True if the dojo_report contains dojo_trial with the given ecut.
        If ecut is None, we test if dojo_trial is present.
        """
        if dojo_trial not in self.ALL_TRIALS:
            raise self.Error("dojo_trial `%s` is not a registered DOJO TRIAL" % dojo_trial)

        if ecut is None:
            return dojo_trial in self
        else:
            #key = self._ecut2key(ecut)
            key = ecut
            try:
                self[dojo_trial][key]
                return True
            except KeyError:
                return False

    def add_ecuts(self, new_ecuts):
        """Add a list of new ecut values."""
        # Be careful with the format here! it should be %.1f
        # Select the list of ecuts reported in the DOJO section.
        prev_ecuts = self["ecuts"]

        for i in range(len(prev_ecuts)-1):
            if prev_ecuts[i] >= prev_ecuts[i+1]:
                raise self.Error("Ecut list is not ordered:\n %s" % prev_ecuts)

        if not isinstance(new_ecuts, Iterable): new_ecuts = [new_ecuts]
        for e in new_ecuts:
            # Find rightmost value less than or equal to x.
            if e < prev_ecuts[0]:
                i = 0
            elif e > prev_ecuts[-1]:
                i = len(prev_ecuts)
            else:
                i = find_le(prev_ecuts, e)
                assert prev_ecuts[i] != e
                i += 1

            prev_ecuts.insert(i, e)

    @property
    def has_hints(self):
        """True if hints on cutoff energy are present."""
        return "hints" in self

    def add_hints(self, hints):
        """Add hints on cutoff energy."""
        hints_dict = {
           "low": {'ecut': hints[0]},
           "normal": {'ecut': hints[1]},
           "high": {'ecut': hints[2]}
                     }
        self["hints"] = hints_dict

    #def validate(self, hints):
    #    Add md5 hash value
    #    self["validated"] = True

    @staticmethod
    def _ecut2key(ecut):
        """Convert ecut to a valid key. ecut can be either a string or a float."""
        if is_string(ecut):
            # Validate string
            i = ecut.index(".")
            if len(ecut[i+1:]) != 1:
                raise ValueError("string %s must have one digit")
            return ecut

        else:
            # Assume float
            return "%.1f" % ecut

    def add_entry(self, dojo_trial, ecut, entry, overwrite=False):
        """
        Add an entry computed with the given ecut to the sub-dictionary associated to dojo_trial.

        Args:
            dojo_trial: String defining the dojo trial.
            ecut: Cutoff energy in Hartree
            entry: Dictionary with data.
            overwrite: By default, this method raises ValueError if this entry is already filled.
        """
        if dojo_trial not in self.ALL_TRIALS:
            raise ValueError("%s is not a registered trial")

        if dojo_trial not in self: self[dojo_trial] = {}
        section = self[dojo_trial]

        key = self._ecut2key(ecut)
        if key in section and not overwrite:
            raise self.Error("Cannot overwrite key %s in dojo_trial %s" % (key, dojo_trial))

        # Add entry to section.
        section[key] = entry

    def find_missing_entries(self):
        """
        Check the DojoReport.
        This function tests if each trial contains an ecut entry.
        Return a dictionary {trial_name: [list_of_missing_ecuts]}
        mapping the name of the Dojo trials to the list of ecut values that are missing
        """
        d = {}

        for trial in self.ALL_TRIALS:
            data = self.get(trial, None)
            if data is None:
                # Gbrv results do not contain noble gases so ignore the error
                if "gbrv" in trial and self.element.is_noble_gas:
                    assert data is None
                    continue
                d[trial] = self.ecuts

            else:
                computed_ecuts = self[trial].keys()
                for e in self.ecuts:
                    if e not in computed_ecuts:
                        if trial not in d: d[trial] = []
                        d[trial].append(e)

        if not d:
            assert len(computed_ecuts) == len(self.ecuts)

        return d

    def get_ecut_dfactprime(self):
        """Return numpy arrays wit ecut list and the corresponding dfactprime values."""
        data = self["deltafactor"]
        ecuts, values= data.keys(), []
        values = np.array([data[e]["dfactprime_meV"] for e in ecuts])
        return np.array(ecuts), values

    def compute_hints(self):
        ecuts, dfacts = self.get_ecut_dfactprime()
        abs_diffs = np.abs((dfacts - dfacts[-1]))
        #print(list(zip(ecuts, dfacts)))
        #print(abs_diffs)

        hints = 3 * [None]
        for ecut, adiff in zip(ecuts, abs_diffs):
            for i in range(3):
                if adiff <= self.ATOLS[i] and hints[i] is None:
                    hints[i] = ecut
                if adiff > self.ATOLS[i]:
                    hints[i] = None
        return hints

    def check(self, check_trials=None):
        """
        Check the dojo report for inconsistencies.
        Return a string with the errors found in the DOJO_REPORT.

        Args:
            check_trials: string or list of strings selecting the trials to be tested.
                If None, all trials are analyzed.
        """
        check_trials = self.ALL_TRIALS if check_trials is None else list_strings(check_trials)
        errors = []
        app = errors.append

        for k in ("version", "ppgen_hints", "md5"):
            if k not in self: app("%s is missing" % k)

        # Check if we have computed each trial for the full set of ecuts in global_ecuts
        global_ecuts = self.ecuts

        missing = defaultdict(list)
        for trial in check_trials:
            for ecut in global_ecuts:
                if not self.has_trial(trial, ecut=ecut):
                    missing[trial].append(ecut)

        if missing:
            app("%s: the following ecut values are missing:" % self.symbol)
            for trial, ecuts in missing.items():
                app("    %s: %s" % (trial, ecuts))

        for trial in check_trials:
            if not self.has_trial(trial): continue
            for ecut in self[trial]:
                if ecut not in global_ecuts:
                    app("%s: ecut %s is not in the global list" % (trial, ecut))

        return "\n".join(errors)

    #def print_table(self, stream=sys.stdout):
    #    from monty.pprint import pprint_table
    #    pprint_table(self.get_dataframe(), out=stream)

    @add_fig_kwargs
    def plot_etotal_vs_ecut(self, ax=None, inv_ecut=False, **kwargs):
        """
        plot the convergence of the total energy as function of the energy cutoff ecut

        Args:
            ax: matplotlib Axes, if ax is None a new figure is created.

        Returns:
            `matplotlib` figure.
        """
        # Extract the total energy of the AE relaxed structure (4).
        d = OrderedDict([(ecut, data["etotals"][4]) for ecut, data in self["deltafactor"].items()])

        # Ecut mesh in Ha
        ecuts = np.array(list(d.keys()))
        ecut_min, ecut_max = np.min(ecuts), np.max(ecuts)

        # Energies per atom in meV and difference wrt 'converged' value
        num_sites = [v["num_sites"] for v in self["deltafactor"].values()][0]
        etotals_mev = np.array([d[e] for e in ecuts]) * 1000  / num_sites
        ediffs = etotals_mev - etotals_mev[-1]

        ax, fig, plt = get_ax_fig_plt(ax)
        #ax.yaxis.set_view_interval(-5, 5)

        lines, legends = [], []

        xs = 1/ecuts if inv_ecut else ecuts
        ys = etotals_mev if inv_ecut else ediffs

        line, = ax.plot(xs, ys, "-o", color="blue") #, linewidth=3.0, markersize=15)
        lines.append(line)

        label = kwargs.pop("label", None)
        if label is not None: ax.legend(lines, [label], loc='best', shadow=True)

        high_hint = self["ppgen_hints"]["high"]["ecut"]
        #ax.vlines(high_hint, min(ediffs), max(ediffs))
        #ax.vlines(high_hint, 0.5, 1.5)
        #ax.scatter([high_hint], [1.0], s=20) #, c='b', marker='o', cmap=None, norm=None)
        #ax.arrow(high_hint, 1, 0, 0.2, head_width=0.05, head_length=0.1, fc='k', ec='k',head_starts_at_zero=False)

        #ax.hlines(5, ecut_min, ecut_max, label="5.0")
        #ax.hlines(1, ecut_min, ecut_max, label="1.0")
        #ax.hlines(0.5, ecut_min, ecut_max, label="0.2")

        # Set xticks and labels.
        ax.grid(True)
        ax.set_xlabel("Ecut [Ha]")
        ax.set_xticks(xs)
        ax.set_ylabel("Delta Etotal/natom [meV]")
        #ax.set_xlim(0, max(xs))

        # Use logscale if possible.
        if all(ediffs[:-1] > 0):
            ax.set_yscale("log")
            ax.set_xlim(xs[0]-1, xs[-2]+1)

        return fig

    @add_fig_kwargs
    def plot_deltafactor_eos(self, ax=None, **kwargs):
        """
        plot the EOS computed with the deltafactor setup.

        Args:
            ax: matplotlib :class:`Axes` or None if a new figure should be created.

        ================  ==============================================================
        kwargs            Meaning
        ================  ==============================================================
        cmap              Color map. default `jet`
        ================  ==============================================================

        Returns:
            `matplotlib` figure.
        """
        ax, fig, plt = get_ax_fig_plt(ax)

        trial = "deltafactor"
        ecuts = self[trial].keys()
        num_ecuts = len(ecuts)

        cmap = kwargs.pop("cmap", None)
        if cmap is None: cmap = plt.get_cmap("jet")

        for i, ecut in enumerate(ecuts):
            d = self[trial][ecut]
            num_sites, volumes, etotals = d["num_sites"], np.array(d["volumes"]), np.array(d["etotals"])

            # Use same fit as the one employed for the deltafactor.
            eos_fit = EOS.DeltaFactor().fit(volumes/num_sites, etotals/num_sites)

            label = "ecut %.1f" % ecut if i % 2 == 0 else ""
            label = "ecut %.1f" % ecut
            eos_fit.plot(ax=ax, text=False, label=label, color=cmap(i/num_ecuts, alpha=1), show=False)

        return fig

    @add_fig_kwargs
    def plot_deltafactor_convergence(self, xc, code="WIEN2k", what=None, ax_list=None, **kwargs):
        """
        plot the convergence of the deltafactor parameters wrt ecut.

        Args:
            xc=String or XcFunc object specifying the XC functional. E.g "PBE" or XcFunc.from_name("PBE"
            code: Reference code
            ax_list: List of matplotlib Axes, if ax_list is None a new figure is created

        Returns:
            `matplotlib` figure.
        """
        all = ["dfact_meV", "dfactprime_meV", "v0", "b0_GPa", "b1"]
        if what is None:
            keys = all
        else:
            what = list_strings(what)
            if what[0].startswith("-"):
                # Exclude keys
                #print([type(w) for w in what])
                what = [w[1:] for w in what]
                keys = [k for k in all if k not in what]
            else:
                keys = what

        # get reference entry
        from pseudo_dojo.refdata.deltafactor import df_database
        reference = df_database(xc=xc).get_entry(symbol=self.symbol, code=code)

        d = self["deltafactor"]
        ecuts = list(d.keys())

        import matplotlib.pyplot as plt
        if ax_list is None:
            fig, ax_list = plt.subplots(nrows=len(keys), ncols=1, sharex=True, squeeze=False)
            ax_list = ax_list.ravel()
        else:
            fig = plt.gcf()

        #ax_list, fig, plt = get_axarray_fig_plt(ax_list, nrows=len(keys), ncols=1, sharex=True, squeeze=False)

        if len(keys) != len(ax_list):
            raise ValueError("len(keys)=%s != len(ax_list)=%s" %  (len(keys), len(ax_list)))

        for i, (ax, key) in enumerate(zip(ax_list, keys)):
            values = np.array([float(d[ecut][key]) for ecut in ecuts])
            #try:
            refval = getattr(reference, key)
            #except AttributeError:
            #    refval = 0.0

            # Plot difference pseudo - ref.
            ax.plot(ecuts, values - refval, "o-")

            ax.grid(True)
            ax.set_ylabel("$\Delta$" + key)
            if i == len(keys) - 1: ax.set_xlabel("Ecut [Ha]")

            if key == "dfactprime_meV":
                # Add horizontal lines (used to find hints for ecut).
                last = values[-1]
                xmin, xmax = min(ecuts), max(ecuts)
                for pad, color in zip(self.ATOLS, ("blue", "red", "violet")):
                    ax.hlines(y=last + pad, xmin=xmin, xmax=xmax, colors=color, linewidth=1, linestyles='dashed')
                    ax.hlines(y=last - pad, xmin=xmin, xmax=xmax, colors=color, linewidth=1, linestyles='dashed')

                # Set proper limits so that we focus on the relevant region.
                ax.set_ylim(last - 1.1*self.ATOLS[0], last + 1.1*self.ATOLS[0])

        return fig

    @add_fig_kwargs
    def plot_gbrv_eos(self, struct_type, ax=None, **kwargs):
        """
        Uses Matplotlib to plot the EOS computed with the GBRV setup

        Args:
            ax: matplotlib :class:`Axes` or None if a new figure should be created.

        ================  ==============================================================
        kwargs            Meaning
        ================  ==============================================================
        cmap              Color map. default `jet`
        ================  ==============================================================

        Returns:
            `matplotlib` figure or None if the GBRV test is not present
        """
        ax, fig, plt = get_ax_fig_plt(ax)

        trial = "gbrv_" + struct_type
        # Handle missing entries: noble gases, Hg ...
        if trial not in self: return None
        ecuts = self[trial].keys()
        num_ecuts = len(ecuts)

        cmap = kwargs.pop("cmap", None)
        if cmap is None: cmap = plt.get_cmap("jet")

        for i, ecut in enumerate(ecuts):
            d = self[trial][ecut]
            volumes, etotals = np.array(d["volumes"]), np.array(d["etotals"])

            eos_fit = EOS.Quadratic().fit(volumes, etotals)
            label = "ecut %.1f" % ecut if i % 2 == 0 else ""
            label = "ecut %.1f" % ecut
            eos_fit.plot(ax=ax, text=False, label=label, color=cmap(i/num_ecuts, alpha=1), show=False)

        return fig

    @add_fig_kwargs
    def plot_gbrv_convergence(self, ax_list=None, **kwargs):
        """
        Uses Matplotlib to plot the convergence of the GBRV parameters wrt ecut.

        Args:
            ax_list: List of matplotlib Axes, if ax_list is None a new figure is created

        Returns:
            `matplotlib` figure.
        """
        import matplotlib.pyplot as plt
        stypes = ("fcc", "bcc")
        if ax_list is None:
            fig, ax_list = plt.subplots(nrows=len(stypes), ncols=1, sharex=True, squeeze=False)
            ax_list = ax_list.ravel()
        else:
            fig = plt.gcf()

        #ax_list, fig, plt = get_axarray_fig_plt(ax_list, nrows=len(stypes), ncols=1, sharex=True, squeeze=False)

        if len(stypes) != len(ax_list):
            raise ValueError("len(stypes)=%s != len(ax_list)=%s" %  (len(stypes), len(ax_list)))

        for i, (ax, stype) in enumerate(zip(ax_list, stypes)):
            trial = "gbrv_" + stype
            d = self[trial]
            ecuts = list(d.keys())
            values = np.array([float(d[ecut]["a0_rel_err"]) for ecut in ecuts])

            ax.grid(True)
            ax.set_ylabel("$\Delta$" + trial + "a0_rel_err")

            # Plot difference pseudo - ref.
            ax.plot(ecuts, values, "bo-")
            #ax.hlines(y=0.0, xmin=min(ecuts), xmax=max(ecuts), color="red")
            if i == len(ax_list) - 1: ax.set_xlabel("Ecut [Ha]")

        return fig

    @add_fig_kwargs
    def plot_phonon_convergence(self, ax_list=None, **kwargs):
        """
        Plot the convergence of the phonon modes wrt ecut.

        Args:
            ax_list: List of matplotlib Axes, if ax_list is None a new figure is created

        Returns:
            `matplotlib` figure.
        """
        d = self["phonon"]
        ecuts = list(d.keys())

        l = [(ecut, float(ecut)) for ecut in ecuts]
        s = sorted(l, key=lambda t: t[1])
        max_ecut = s[-1][0]
        s_ecuts = [ecut[0] for ecut in s]

        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(nrows=2, sharex=True)
        #ax_list, fig, plt = get_axarray_fig_plt(ax_list, nrows=len(keys), ncols=1, sharex=True, squeeze=False)

        fmin, fmax = np.inf, -np.inf
        for i, v in enumerate(d[ecuts[0]]):
            values1 = np.array([float(d[ecut][i]) for ecut in s_ecuts])
            fmin = min(fmin, values1.min())
            fmax = max(fmax, values1.max())

            ax[0].plot(s_ecuts, values1, "o-")
            ax[0].grid(True)
            ax[0].set_ylabel("phonon modes [meV] (asr==2)")
            ax[0].set_xlabel("Ecut [Ha]")

            values2 = np.array([float(d[ecut][i]) - float(d[max_ecut][i]) for ecut in s_ecuts])

            ax[1].plot(s_ecuts, values2, "o-")
            ax[1].grid(True)
            ax[1].set_ylabel("w - w(ecut_max) [meV]")
            ax[1].set_xlabel("Ecut [Ha]")

        # Adjust limits.
        fmin -= 10
        fmax += 10
        ax[0].set_ylim(fmin, fmax)

        return fig

    @add_fig_kwargs
    def plot_ebands(self, ecut=None, **kwargs):
        if ecut is None:
            ecut = self['ebands'].keys()[0]
        path = self['ebands'][self['ebands'].keys()[0]]['GSR-nc']

        from abipy.abilab import abiopen
        with abiopen(path) as gsr:
            ebands = gsr.ebands
            fig = ebands.plot_with_edos(ebands.get_edos(width=0.05, step=0.02))
            return fig


from pandas import DataFrame

class DojoDataFrame(DataFrame):
    """Extends pandas DataFrame adding helper functions."""

    ALL_ACCURACIES = ("low", "normal", "high")

    ALL_TRIALS = (
        "ecut",
        "deltafactor",
        "gbrv_bcc",
        "gbrv_fcc",
        "phonon",
        "phwoa",
        "ebands"
    )

    _TRIALS2KEY = {
        "ecut": "ecut",
        "deltafactor": "dfact_meV",
        "gbrv_bcc": "gbrv_bcc_a0_rel_err",
        "gbrv_fcc": "gbrv_fcc_a0_rel_err",
        "phonon": "all",
        "phwoa": "all",
        "ebands": "all"
    }

    _TRIALS2YLABEL = {
        "ecut": "Ecut [Ha]",
        "deltafactor": "$\Delta$-factor [meV]",
        "gbrv_bcc": "BCC $\Delta a_0$ (%)",
        "gbrv_fcc": "FCC $\Delta a_0$ (%)",
        "phonon": "Phonons with ASR",
        "phwoa": "Phonons without ASR",
        "ebands": "Electronic band structure"
    }

    ACC2PLTOPTS = dict(
        low=dict(color="red"),
        normal=dict(color="blue"),
        high=dict(color="green"),
    )

    for v in ACC2PLTOPTS.values():
        v.update(linewidth=2, linestyle='dashed', marker='o', markersize=8)

    def tabulate(self, columns=None, stream=sys.stdout):
        if columns is None:
            accuracies = self.ALL_ACCURACIES
            columns = [acc + "_dfact_meV" for acc in accuracies]
            columns += [acc + "_ecut" for acc in accuracies]
            columns += [acc + "_gbrv_fcc_a0_rel_err" for acc in accuracies]
            columns += [acc + "_gbrv_bcc_a0_rel_err" for acc in accuracies]

        #return self[columns].to_html()
        tablefmt = "grid"
        floatfmt=".2f"
        stream.write(tabulate(self[columns], headers="keys", tablefmt=tablefmt, floatfmt=floatfmt))

    def get_accuracy(self, accuracy):
        columns = [c for c in self if c.startswith(accuracy)]
        return self.__class__(data=self[columns])

    def get_trials(self, accuracies="all"):
        accuracies = self.ALL_ACCURACIES if accuracies == "all" else list_strings(accuracies)

        columns = [acc + "_dfact_meV" for acc in accuracies]
        columns += [acc + "_ecut" for acc in accuracies]
        columns += [acc + "_gbrv_fcc_a0_rel_err" for acc in accuracies]
        columns += [acc + "_gbrv_bcc_a0_rel_err" for acc in accuracies]
        return self.__class__(data=self[columns])

    def select_rows(self, rows):
        if not isinstance(rows, (list, tuple)): rows = [rows]

        data = []
        for index, entry in self.iterrows():
            element = Element.from_Z(entry.Z)
            if element.row in rows:
                data.append(entry)

        return self.__class__(data=data)

    def select_family(self, family):
        data = []
        for index, entry in self.iterrows():
            element = Element.from_Z(entry.Z)
            # e.g element.is_alkaline
            if getattr(element, "is_" + family):
                data.append(entry)
        return self.__class__(data=data)

    @add_fig_kwargs
    def plot_hist(self, what="dfact_meV", bins=400, **kwargs):
        import matplotlib.pyplot as plt
        fig, ax_list = plt.subplots(nrows=len(self.ALL_ACCURACIES), ncols=1, sharex=True, sharey=False, squeeze=True)

        for acc, ax in zip(self.ALL_ACCURACIES, ax_list):
            col = acc + "_" + what
            #print(col)
            #self[col].hist(ax=ax, bins=bins, label=col)
            self[col].plot(ax=ax, kind="bar", label=col)

        return fig

    @add_fig_kwargs
    def plot_trials(self, trials="all", accuracies="all", **kwargs):
        import matplotlib.pyplot as plt
        trials = self.ALL_TRIALS if trials == "all" else list_strings(trials)
        accuracies = self.ALL_ACCURACIES if accuracies == "all" else list_strings(accuracies)

        fig, ax_list = plt.subplots(nrows=len(trials), ncols=1, sharex=True, sharey=False, squeeze=True)

        # See also http://matplotlib.org/examples/pylab_examples/barchart_demo.html
        for i, (trial, ax) in enumerate(zip(trials, ax_list)):
            what = self._TRIALS2KEY[trial]
            ax.set_ylabel(self._TRIALS2YLABEL[trial])
            minval, maxval = np.inf, -np.inf
            for acc in accuracies:
                col = acc + "_" + what
                legend = i == 0
                data = self[col]
                minval, maxval = min(minval, data.min()), max(maxval, data.max())
                data.plot(ax=ax, legend=legend, use_index=True, label=acc, **self.ACC2PLTOPTS[acc])
                #data.plot(ax=ax, kind="bar")

                if i == 0:
                    ax.legend(loc='best', shadow=True, frameon=True) #fancybox=True)

            ax.set_xticks(range(len(data.index)))
            ax.set_xticklabels(data.index)
            #ax.set_xticklabels([root for root, ext in map(os.path.splitext, data.index)])

            # Set ylimits
            #stepsize = None
            #if "gbrv" in trial:
            #    ax.hlines(0.0, 0, len(data.index))
            #    #start, end = -0.6, +0.6
            #    start, end = max(-0.6, minval), min(+0.6, maxval)
            #    if end - start < 0.05: end = start + 0.1
            #    ax.set_ylim(start, end)
            #    ax.yaxis.set_ticks(np.arange(start, end, 0.05))

            if trial == "deltafactor":
                #start, end = 0.0, 15
                start, end  = 0.0, min(15, maxval)
                ax.set_ylim(start, end)
                #ax.yaxis.set_ticks(np.arange(start, end, 0.1))

            #if stepsize is not None:
            #    start, end = ax.get_ylim()
            #    ax.yaxis.set_ticks(np.arange(start, end, stepsize))

            plt.setp(ax.xaxis.get_majorticklabels(), rotation=25)

        return fig
