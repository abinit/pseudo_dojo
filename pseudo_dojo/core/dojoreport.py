"""
"""
import sys
import os
import json
import logging
import copy
import numpy as np
import pandas as pd

from collections import OrderedDict, defaultdict, Iterable
from tabulate import tabulate
from monty.json import MontyEncoder
from monty.string import list_strings, is_string
from monty.termcolor import cprint
from monty.bisect import find_le
from pymatgen.core.periodic_table import Element
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt
from pseudo_dojo.refdata.deltafactor import df_database, df_compute
from pseudo_dojo.refdata.gbrv import gbrv_database
from pseudo_dojo.refdata.lantanides.database import raren_database
from pseudo_dojo.util.dojo_eos import EOS


logger = logging.getLogger(__name__)


def dojo_dfact_results(pseudo, num_sites, volumes, etotals):
    """
    This function computes the deltafactor and returns the dictionary to be inserted
    in the dojoreport file.

    Args:
        pseudo: Pseudopotential object.
        num_sites: Number of sites in unit cell
        volumes: List with unit cell volumes in Ang**3
        etotals: List of total energies in eV.

    Return:
        (dojo_entry, eos_fit)
        where dojo_entry is the Dictionary with results to be inserted in the djrepo file.
        eos_fit is the object storing the results of the EOS fit.
    """
    nan = float('NaN')

    dojo_entry = dict(
        etotals=list(etotals),
        volumes=list(volumes),
        num_sites=num_sites,
        dfact_meV=nan,
        dfactprime_meV=nan,
        v0=nan,
        b0=nan,
        b0_GPa=nan,
        b1=nan,
    )

    volumes, etotals = np.asarray(volumes), np.asarray(etotals)
    eos_fit = None
    try:
        # Use same fit as the one employed for the deltafactor.
        eos_fit = EOS.DeltaFactor().fit(volumes/num_sites, etotals/num_sites)

        # Get reference results (Wien2K).
        wien2k = df_database(pseudo.xc).get_entry(pseudo.symbol)

        # Compute deltafactor estimator.
        dfact = df_compute(wien2k.v0, wien2k.b0_GPa, wien2k.b1,
                           eos_fit.v0, eos_fit.b0_GPa, eos_fit.b1, b0_GPa=True)

        dfactprime_meV = dfact * (30 * 100) / (eos_fit.v0 * eos_fit.b0_GPa)

        dfres = {
            "dfact_meV": dfact,
            "dfactprime_meV": dfactprime_meV,
            "v0": eos_fit.v0,
            "b0": eos_fit.b0,
            "b0_GPa": eos_fit.b0_GPa,
            "b1": eos_fit.b1,
        }

        for k, v in dfres.items():
            v = v if not isinstance(v, complex) else nan
            dfres[k] = v

        dojo_entry.update(dfres)

    except EOS.Error as exc:
        dojo_entry["_exceptions"] = str(exc)

    return dojo_entry, eos_fit


def dojo_gbrv_results(pseudo, struct_type, num_sites, volumes, etotals):
    """
    This function computes the GBRV results and returns the dictionary
    to be inserted in the dojoreport file.

    Args:
        pseudo: Pseudopotential object.
        struct_type: "fcc" or "bcc"
        num_sites: Number of sites in unit cell
        volumes: List with unit cell volumes in Ang**3
        etotals: List of total energies in eV.

    Return:
        (dojo_entry, eos_fit)
        where dojo_entry is the Dictionary with results to be inserted in the djrepo file.
        eos_fit is the object storing the results of the EOS fit.
    """
    # Read etotals and fit E(V) with a parabola to find the minimum
    assert len(etotals) == len(volumes)

    dojo_entry = dict(
        volumes=list(volumes),
        etotals=list(etotals),
        num_sites=num_sites,
    )

    eos_fit = None
    try:
        eos_fit = EOS.Quadratic().fit(volumes, etotals)
    except EOS.Error as exc:
        dojo_entry["_exceptions"] = str(exc)
        return dojo_entry, eos_fit

    # Function to compute cubic a0 from primitive v0 (depends on struct_type)
    vol2a = {"fcc": lambda vol: (4 * vol) ** (1/3.),
             "bcc": lambda vol: (2 * vol) ** (1/3.),
             }[struct_type]

    a0 = vol2a(eos_fit.v0)
    dojo_entry.update(dict(
        v0=eos_fit.v0,
        b0=eos_fit.b0,
        #b1=eos_fit.b1, # infinity
        a0=a0,
        struct_type=struct_type
    ))

    db = gbrv_database(pseudo.xc)
    ref = db.get_entry(pseudo.symbol, stype=struct_type)

    pawabs_err = a0 - ref.gbrv_paw
    pawrel_err = 100 * (a0 - ref.gbrv_paw) / ref.gbrv_paw

    # AE results for P and Hg are missing.
    if ref.ae is not None:
        abs_err = a0 - ref.ae
        rel_err = 100 * (a0 - ref.ae) / ref.ae
    else:
        # Use GBRV_PAW as reference.
        abs_err = pawabs_err
        rel_err = pawrel_err

    print("for GBRV struct_type: ", struct_type, "a0= ", a0, "Angstrom")
    print("AE - THIS: abs_err = %f, rel_err = %f %%" % (abs_err, rel_err))
    print("GBRV-PAW - THIS: abs_err = %f, rel_err = %f %%" % (pawabs_err, pawrel_err))

    dojo_entry["a0_abs_err"] = abs_err
    dojo_entry["a0_rel_err"] = rel_err

    return dojo_entry, eos_fit


class DojoReportError(Exception):
    """Exception raised by DoJoReport."""


class DojoReport(dict):
    """
    Dict-like object with the validation results.
    This object is usually created via the class methods:

        DojoReport.from_file and DojoReport.empty_from_pseudo.

    {
    "version": "1.0"
    "symbol": "H",
    "pseudo_type": "NC",
    "md5": "13198abb7506a840b7d46ef46b54d789",
    "ppgen_hints": {
        "low": {"ecut": 30.0,  "pawecutdg": 30.0},
        "normal": {"ecut": 34.0, "pawecutdg": 34.0},
        "high": {"ecut": 39.0, "pawecutdg": 39.0}
    },
    "hints": {
        "low": {"ecut": 30.0,  "pawecutdg": 30.0},
        "normal": {"ecut": 34.0, "pawecutdg": 34.0},
        "high": {"ecut": 39.0, "pawecutdg": 39.0},
    },
    "ecuts": [29.0, 31.0, 33.0],
    "deltafactor": {}
    "gbrv_bcc": {},
    "gbrv_fcc": {},
    "ghosts": []
    "phonons": []
    }

    "ecut": 32.0
    "pawecutdg": 64.0,
    "b0": 0.06400805819081799,
    "b0_GPa": 10.255221080448488,
    "b1": 2.6449207740813594,
    "dfact_meV": 0.2774768889565598,
    "dfactprime_meV": 4.701668998922405,
    "etotals": []
    "num_sites": 4,
    "v0": 17.264380250637252,
    "volumes": [],
    """

    # List of dojo_trials. Remember to update the list if you add a new test to the DOJO_REPORT
    ALL_TRIALS = [
        "deltafactor",
        "gbrv_bcc",
        "gbrv_fcc",
        "phgamma",
        "ghosts",
    ]

    # Add trials done with SOC.
    ALL_TRIALS += [n + "_soc" for n in ALL_TRIALS]

    _TRIALS2KEY = {
        "deltafactor": "dfact_meV",
        "gbrv_bcc": "a0_rel_err",
        "gbrv_fcc": "a0_rel_err",
        "phgamma": "all",
        #"ghosts": "all",
    }

    # We use three different level of accuracy.
    ALL_ACCURACIES = ("low", "normal", "high")
    ACC2COLOR = {"low": "yellow", "normal": "green", "high": "red"}

    # Tolerances on the deltafactor prime (in eV) used for the hints.
    #ATOLS = (0.5, 0.1, 0.02)
    ATOLS = (0.5, 0.3, 0.1)
    # For noble gasses:
    #ATOLS = (1.0, 0.2, 0.04)

    # Version of the DojoReport.
    LAST_VERSION = "1.0"

    Error = DojoReportError

    @classmethod
    def from_file(cls, filepath):
        """Read the DojoReport from file."""
        filepath = os.path.abspath(filepath)
        with open(filepath, "rt") as fh:
            d = json.load(fh)
            new = cls(**d)
            new.path = filepath
            # TODO
            #new.xc = XcFunc.from_dict(new["xc"])
            return new

    @classmethod
    def empty_from_pseudo(cls, pseudo, ppgen_hints, devel=False):
        """
        Initialize an empty `DojoReport` from the pseudo and an initial guess
        for the cutoff energies in Hartree

        Args:
            pseudo: :class:`Pseudo` object.
            ppgen_hints: Initial hints on the cutoff energy provided by the pp generator.
                Dictionary [accuracy]["ecut"] --> ecut_value
        """
        # Build initial list of cutoff energies for tests.
        #dense_right = np.arange(ppgen_ecut, ppgen_ecut + 6*2, step=2)
        #dense_left = np.arange(max(ppgen_ecut-6, 2), ppgen_ecut, step=2)
        #coarse_high = np.arange(ppgen_ecut + 15, ppgen_ecut + 35, step=5)

        new = cls()

        estart = ppgen_hints["high"]["ecut"]
        dense_right = np.linspace(estart - 10, estart + 10, num=11)
        ecuts = list(dense_right) + [dense_right[-1] + 8, dense_right[-1] + 10,]

        # devel is for tuning the pseudo, only two cutoffs
        # development run: few, relatively high ecut calculations
        if devel: ecuts = [estart, estart + 2]

        if pseudo.isnc:
            pseudo_type = "NC"
        elif pseudo.ispaw:
            pseudo_type = "PAW"
        else:
            raise TypeError("Neither NC nor PAW pseudo!")

        new.update(
            basename=pseudo.basename,
            version=cls.LAST_VERSION,
            symbol=pseudo.symbol,
            pseudo_type=pseudo_type,
            xc=pseudo.xc.as_dict(),
            md5=pseudo.compute_md5(),
            ppgen_hints=ppgen_hints,
            ecuts=ecuts,
        )
        new.path = pseudo.djrepo_path

        return new

    # TODO Remove
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

        # TODO
        #new.xc = XcFunc.from_dict(new["xc"])

    # TODO Remove
    def reorder(self):
        for trial in self.ALL_TRIALS:
            # Convert ecut to float and build an OrderedDict (results are indexed by ecut in ascending order)
            try:
                d = self[trial]
                #print("Reordering", trial)
            except KeyError:
                continue
            ecuts_keys = sorted([(float(k), k) for k in d], key=lambda t: t[0])
            self[trial] = OrderedDict([(t[0], d[t[1]]) for t in ecuts_keys])

    def __str__(self):
        """String representation."""
        return(json.dumps(self, indent=4))

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
        """list of ecuts that should be present in the dojo_trial sub-dicts"""
        return self["ecuts"]

    @property
    def trials(self):
        """Set of strings with the trials present in the report."""
        return set(list(self.keys())).intersection(self.ALL_TRIALS)

    @property
    def exceptions(self):
        """List of exceptions."""
        return self.get("_exceptions", [])

    def push_exception(self, string):
        """Register an exception."""
        if "_exceptions" not in self:
            self["_exceptions"] = []
        self["_exceptions"].append(string)

    def remove_exceptions(self):
        """Remove the exception entry from the dictionary."""
        return self["_exceptions"].pop()

    def deepcopy(self):
        """Deepcopy of the object."""
        return copy.deepcopy(self)

    def json_write(self, filepath=None):
        """
        Write data to file in JSON format.
        If filepath is None, self.path is used.
        """
        filepath = self.path if filepath is None else filepath
        with open(filepath, "wt") as fh:
            #json.dump(self, fh, indent=-1, sort_keys=True)
            json.dump(self, fh, indent=-1, sort_keys=True, cls=MontyEncoder)

    def has_trial(self, dojo_trial, ecut=None):
        """
        True if the dojo_report contains dojo_trial with the given ecut.
        If ecut is None, we test if dojo_trial is present.
        """
        if dojo_trial not in self: return False
        if ecut is None: return dojo_trial in self

        # input ecut could be either float or string
        if ecut in self[dojo_trial]: return True
        ecut_str = self._ecut2key(ecut)
        if ecut_str in self[dojo_trial]: return True
        return False

    def remove_trial(self, dojo_trial, ecut=None, write=False):
        """
        Remove the entry associated to `dojo_trial` and write new JSON file.
        If ecut is None, the entire `dojo_trial` is removed.
        The writing of the JSON file can be postponed by setting write=False.
        """
        if ecut is None:
            res = self.pop(dojo_trial)
        else:
            if not is_string(ecut): ecut = "%.1f" % ecut
            res = self[dojo_trial].pop(ecut)

        if write: self.json_write()
        return res

    def get_pdframe(self, dojo_trial, *args):
        """
        Return a pandas `DataFrame` with the results for this particular dojo_trial.
        The frame is sorted according to the value of ecut.

        This function is used to plot/analyze data for a single pseudo.
        Args is the list of variable that should be used to construct the frame a.k.a column names.
        If args is empty, all the keys found in the dictionary are used

        Examples:

        To extract all data available:

            frame = dojo_report.get_pdframe("deltafactor")

        To build a frame with ["ecut", "dfact_meV", "b0"]:

            frame = dojo_report.get_pdframe("deltafactor", "dfact_meV", "b0")
        """
        if dojo_trial not in self:
            raise self.Error("dojo_trial %s not present" % dojo_trial)

        if not args:
            #args = list(self[dojo_trial][0].keys())
            #OLD DOJO
            args = []
            for ecut, data in self[dojo_trial].items():
                args.extend(list(data.keys()))
            args = set(args)

        # Build list of dictionaries.
        dict_list = []

        # This code is compatible with dict: ecut --> entry
        for ecut, data in self[dojo_trial].items():
            d = dict(ecut=float(ecut), pawecutdg=data.get("pawecutdg", None))
            d.update({k: data[k] for k in args})
            dict_list.append(d)

        # TODO: This code assumes a list of entries instead of a dict: ecut --> entry
        #for entry in self[dojo_trial]:
        #    d = entry["ecut"]
        #    d["pawecutdg"] = entry.get("pawecutdg", None)
        #    d.update({k: entry[k] for k in args})
        #    dict_list.append(d)

        # Build DataFrame from dict_list and sort wrt ecut.
        frame = pd.DataFrame(dict_list)
        return frame.sort_values("ecut")

    def add_ecuts(self, new_ecuts, write=False):
        """
        Add a list of new set of ecuts to the global list
        If write is True, the JSON file is immediately updated with the new data.
        """
        if not isinstance(new_ecuts, Iterable): new_ecuts = [new_ecuts]

        # Be careful with the format here! it should be %.1f
        # Select the list of ecuts reported in the DOJO section.
        prev_ecuts = self["ecuts"]

        # prev_ecuts should be ordered.
        for i in range(len(prev_ecuts)-1):
            if prev_ecuts[i] >= prev_ecuts[i+1]:
                raise self.Error("Ecut list is not ordered:\n %s" % prev_ecuts)

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

        if write: self.json_write()

    @property
    def has_hints(self):
        """True if hints on the cutoff energy are present."""
        return "hints" in self

    #@property
    #def hints(self):

    def add_hints(self, hints, write=False):
        """
        Add hints on the cutoff energy.
        If write is True, the JSON file is immediately updated with the new data.
        """
        self["hints"] = {
           "low": {'ecut': hints[0]},
           "normal": {'ecut': hints[1]},
           "high": {'ecut': hints[2]}
        }

        if write: self.json_write()

    @property
    def isvalidated(self):
        """True if the dojoreport has been validated."""
        return "validation" in self and self.has_hints

    def ipw_validate(self):
        """
        Return an ipython widget with controllers to validate the pseudo.
        """
        import ipywidgets as ipw
        low_ecut = ipw.FloatText(description='Low ecut:')
        normal_ecut = ipw.FloatText(description='Normal ecut:')
        high_ecut = ipw.FloatText(description='High ecut:')
        new_validation = ipw.Checkbox(description="New validation", value=False)
        s = " "
        if "tags" in self: s = ",".join(self["tags"])
        tags = ipw.Text(description='Tags:', value=s)
        #validated_by = ipw.Text(description="Validated by:")
        #validated_by.value = "authors"
        ok_button = ipw.Button(description="Validate")
        if self.has_hints:
            low_ecut.value = self["hints"]["low"]["ecut"]
            normal_ecut.value = self["hints"]["normal"]["ecut"]
            high_ecut.value = self["hints"]["high"]["ecut"]

        def on_button_clicked(b):
            """Callback called to validate the dojo report."""
            print(low_ecut.value, normal_ecut.value, high_ecut.value)
            if not low_ecut.value <= normal_ecut.value <= high_ecut.value:
                raise ValueError("not low_ecut.value <= normal_ecut.value <= high_ecut.value")
            #if not validated_by.value:
            #    raise ValueError("validated_by field must be filled")
            if self.isvalidated and not new_validation.value:
                raise ValueError("DojoReport is already validated. Use new_validation")

            from time import gmtime, strftime
            self['validation'] = {
                #'validated_by': validated_by.value,
                'validated_on': strftime("%Y-%m-%d %H:%M:%S", gmtime())
            }

            # Add hints
            hints = {"low": {}, "normal": {}, "high": {}}
            hints["low"]["ecut"] = low_ecut.value
            hints["normal"]["ecut"] = normal_ecut.value
            hints["high"]["ecut"] = high_ecut.value
            self["hints"] = hints
            self["tags"] = tags.value.split(",")
            print(hints, tags.value.split(","))
            self.json_write()

        ok_button.on_click(on_button_clicked)
        return ipw.VBox(children=[low_ecut, normal_ecut, high_ecut, new_validation, tags, ok_button])

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

    def add_entry(self, dojo_trial, ecut, entry, overwrite=False, write=False):
        """
        Add an entry computed with the given ecut to the sub-dictionary associated to dojo_trial.

        Args:
            dojo_trial: String defining the dojo trial.
            ecut: Cutoff energy in Hartree
            entry: Dictionary with data.
            overwrite: By default, this method raises ValueError if this entry is already filled.
            write: if True, the JSON file is immediately updated with the new data.
        """
        if dojo_trial not in self.ALL_TRIALS and dojo_trial.replace("_soc", "", 1) not in self.ALL_TRIALS:
            raise ValueError("%s is not a registered trial")

        if dojo_trial not in self: self[dojo_trial] = {}
        section = self[dojo_trial]

        key = self._ecut2key(ecut)
        if key in section and not overwrite:
            raise self.Error("Cannot overwrite key %s in dojo_trial %s" % (key, dojo_trial))

        # Add entry to section.
        section[key] = entry

        if write: self.json_write()

    def find_missing_entries(self):
        """
        This function tests if each trial contains an ecut entry.
        Return a dictionary {trial_name: [list_of_missing_ecuts]}
        mapping the name of the Dojo trials to the list of ecut values that are missing
        """
        misse = {}

        if "ghosts" not in self: misse["ghosts"] = [0]

        for trial in self.ALL_TRIALS:
            data = self.get(trial, None)
            if data is None:
                # Gbrv results do not contain noble gases so ignore the error
                if "gbrv" in trial and self.element.is_noble_gas:
                    assert data is None
                    continue
                misse[trial] = self.ecuts

            else:
                computed_ecuts = list(data.keys())
                for e in self.ecuts:
                    if e not in computed_ecuts:
                        if trial not in misse: misse[trial] = []
                        misse[trial].append(e)

        if not misse:
            assert len(computed_ecuts) == len(self.ecuts)

        return misse

    def get_ecut_dfactprime(self):
        """Return numpy arrays wit ecut list and the corresponding dfactprime values."""
        data = self["deltafactor"]
        ecuts, values = data.keys(), []
        values = np.array([data[e]["dfactprime_meV"] for e in ecuts])
        return np.array(ecuts), values

    def get_last_df_results(self, with_soc=False):
        """
        Return dictionary with the last value i.e. the best estimate of deltafactor
        and deltafactor_prime. Empty dictionary if results are not available

        Args:
            with_soc: If True, the results obtained with SOC are returned (if available).
            In this case, the name of variable contains the `_soc` suffix at the end e.g.
            `v0` becomes `v0_soc`.
        """
        trial = "deltafactor" if not with_soc else "deltafactor_soc"
        try:
            data = self[trial]
        except KeyError:
            return {}

        # Get the values associated with the last ecut (highest value).
        vnames = ["dfact_meV", "dfactprime_meV", "v0", "b0_GPa", "b1"]
        frame = self.get_pdframe(trial, *vnames)
        d = {vname: frame[vname].iloc[-1] for vname in vnames}

        # Add `_soc` prefix.
        if with_soc:
            d = {k + "_soc": d[k] for k in d}

        return d

    def get_last_gbrv_results(self, struct_type, with_soc=False):
        """
        Return dictionary with the last value i.e. the best estimate of the
        gbrv results. Empty dictionary if results are not available

        Args:
            struct_type: "bcc" or "fcc".
            with_soc: If True, the results obtained with SOC are returned (if available).
        """
        trial = "gbrv_" + struct_type
        if with_soc: trial += "_soc"
        try:
            data = self[trial]
        except KeyError:
            return {}

        # Get the values associated with the last ecut (highest value).
        vnames = ["a0", "a0_abs_err", "a0_rel_err"]
        frame = self.get_pdframe(trial, *vnames)
        d = {trial + "_" + vname: frame[vname].iloc[-1] for vname in vnames}

        # Add `_soc` prefix.
        #if with_soc:
        #    d = {k + "_soc": d[k] for k in d}

        return d

    def check(self, check_trials=None):
        """
        Check the dojo report for inconsistencies.
        Return a string with the errors found in the DOJO_REPORT.

        Args:
            check_trials: string or list of strings selecting the trials to be tested.
                If None, all trials are analyzed.
        """
        check_trials = self.ALL_TRIALS if check_trials is None else list_strings(check_trials)
        # Relativistic pseudos (_r.psp8) must have trials done with/without SOC
        if "_r" not in self["basename"]:
            check_trials = [t for t in check_trials if not t.endswith("_soc")]

        errors = []
        app = errors.append
        for k in ("version", "ppgen_hints", "md5", "ghosts"):
            if k not in self: app("%s is missing" % k)

        # Check if we have computed each trial for the full set of ecuts in global_ecuts
        global_ecuts = self.ecuts

        # TODO: report should contain XC
        missing = defaultdict(list)
        for trial in check_trials:
            # GBRV results do not contain noble gases, Hg and Po so ignore the error
            if "gbrv" in trial and (self.element.is_noble_gas or
                                    self.element.is_lanthanoid or
                                    self.element.is_actinoid or
                                    self.symbol in ("Hg", "Po")): continue

            if "deltafactor" in trial and self.symbol != "Lu" and \
                (self.element.is_lanthanoid or self.element.is_actinoid): continue

            if trial not in self:
                missing[trial].append(0)

            #for ecut in global_ecuts:
            #    if not self.has_trial(trial, ecut=ecut):
            #        missing[trial].append(ecut)

        if missing:
            app("%s: the following ecut values are missing:" % self.symbol)
            for trial, ecuts in missing.items():
                app("    %s: %s" % (trial, ecuts))

        #for trial in check_trials:
        #    if not self.has_trial(trial): continue
        #    for ecut in self[trial]:
        #        if ecut not in global_ecuts:
        #            app("%s: ecut %s is not in the global list" % (trial, ecut))

        if "hints" not in self:
            app("hints are missing")
        else:
            for acc in ["low", "normal", "high"]:
                assert self["hints"][acc]["ecut"] > 0

        for k in ("md5", "md5_upf", "md5_psml"):
            if k not in self:
                app("%s not djrepo" % k)

        return "\n".join(errors)

    ##################
    # Plotting tools #
    ##################

    @add_fig_kwargs
    def plot_etotal_vs_ecut(self, ax=None, inv_ecut=False, with_soc=False, **kwargs):
        """
        Plot the convergence of the total energy as function of the energy cutoff ecut

        Args:
            ax: matplotlib Axes, if ax is None a new figure is created.
            inv_ecut: True to plot etotal vs 1/ecut.
            with_soc: If True, the results obtained with SOC are plotted (if available).

        Returns:
            `matplotlib` figure or None if the deltafactor test is not present
        """
        trial = "deltafactor" if not with_soc else "deltafactor_soc"
        if trial not in self:
            cprint("dojo report does not contain trial: %s" % str(trial), "red")
            return None

        label = kwargs.pop("label", None)

        # Get Data Frame
        frame = self.get_pdframe(trial, "etotals", "num_sites")

        # Get Ecut mesh in Ha
        ecuts = np.array(frame["ecut"])
        ecut_min, ecut_max = np.min(ecuts), np.max(ecuts)

        # Extract the total energy of the AE relaxed structure (4).
        minenes = np.array([etots[4] for etots in frame["etotals"]])
        num_sites = np.array(frame["num_sites"])
        assert np.all(num_sites == num_sites[0])
        num_sites = num_sites[0]

        # Energies per atom in meV and difference wrt 'converged' value
        etotals_mev = minenes * 1000 / num_sites
        ediffs = etotals_mev - etotals_mev[-1]

        ax, fig, plt = get_ax_fig_plt(ax)
        #ax.yaxis.set_view_interval(-5, 5)

        lines, legends = [], []
        xs = 1/ecuts if inv_ecut else ecuts
        ys = etotals_mev if inv_ecut else ediffs

        line, = ax.plot(xs, ys, "-o", color="blue") #, linewidth=3.0, markersize=15)
        lines.append(line)

        # Add vertical lines at hints.
        if self.has_hints:
            vmax = ys.max()
            vmin = ys.min() if inv_ecut else np.min([y for y in ys if y > 0])
            for acc in self.ALL_ACCURACIES:
                x0 = self["hints"][acc]["ecut"]
                if inv_ecut: x0 = 1/x0
                ax.vlines(x0, vmin+0.001, vmax, colors=self.ACC2COLOR[acc], linestyles="dashed")

        if label is not None: ax.legend(lines, [label], loc='best', shadow=True)

        # Set xticks and labels.
        ax.grid(True)
        ax.set_xticks(xs)
        if not inv_ecut:
            ax.set_xlabel("Ecut [Ha]")
            ax.set_ylabel("Delta Etotal/natom [meV]")
        else:
            ax.set_xlabel("1/Ecut [Ha]")
            ax.set_ylabel("Etotal/natom [meV]")
            ax.set_ylim(ys.min() - 0.5, ys.min() + 2.5)
            if self.has_hints:
                x0 = self["hints"]["low"]["ecut"] - 4
                ax.set_xlim(0, 1/x0)

        # Use logscale if possible.
        if not inv_ecut and all(ediffs[:-1] > 0):
            ax.set_yscale("log")
            ax.set_xlim(xs[0]-1, xs[-2]+1)

        plt.tight_layout()

        return fig

    @add_fig_kwargs
    def plot_deltafactor_eos(self, ax=None, with_soc=False, **kwargs):
        """
        Plot the equation of state computed with the deltafactor setup.

        Args:
            ax: matplotlib :class:`Axes` or None if a new figure should be created.
            with_soc: If True, the results obtained with SOC are plotted (if available).

        ================  ==============================================================
        kwargs            Meaning
        ================  ==============================================================
        cmap              Color map. default `jet`
        ================  ==============================================================

        Returns:
            `matplotlib` figure or None if the deltafactor test is not present
        """
        trial = "deltafactor" if not with_soc else "deltafactor_soc"
        if trial not in self:
            cprint("dojo report does not contain trial: %s" % str(trial), "red")
            return None

        ax, fig, plt = get_ax_fig_plt(ax)
        cmap = kwargs.pop("cmap", plt.get_cmap("jet"))

        # Get DataFrame.
        frame = self.get_pdframe(trial, "num_sites", "volumes", "etotals")
        ecuts = frame["ecut"]
        num_sites = np.array(frame["num_sites"])
        assert np.all(num_sites == num_sites[0])
        num_sites = num_sites[0]

        for i, ecut in enumerate(ecuts):
            # Subframe with this value of ecut.
            ecut_frame = frame.loc[frame["ecut"] == ecut]
            assert ecut_frame.shape[0] == 1
            # Extract volumes and energies for this ecut.
            volumes = (np.array(list(ecut_frame["volumes"].values), dtype=np.float)).flatten()
            etotals = (np.array(list(ecut_frame["etotals"].values), dtype=np.float)).flatten()

            # Use same fit as the one employed for the deltafactor.
            eos_fit = EOS.DeltaFactor().fit(volumes/num_sites, etotals/num_sites)
            eos_fit.plot(ax=ax, text=False, label="ecut %.1f" % ecut, color=cmap(i/len(ecuts), alpha=1), show=False)

        plt.tight_layout()

        return fig

    @add_fig_kwargs
    def plot_deltafactor_convergence(self, xc, code="WIEN2k", with_soc=False, what=None, ax_list=None, **kwargs):
        """
        Plot the convergence of the deltafactor parameters wrt ecut.

        Args:
            xc: String or XcFunc object specifying the XC functional. E.g "PBE" or XcFunc.from_name("PBE")
            code: Reference code.
            with_soc: If True, the results obtained with SOC are plotted (if available).
            what:
            ax_list: List of matplotlib Axes, if ax_list is None a new figure is created

        Returns:
            `matplotlib` figure or None if the deltafactor test is not present
        """
        trial = "deltafactor" if not with_soc else "deltafactor_soc"
        if trial not in self:
            cprint("dojo report does not contain trial: %s" % str(trial), "red")
            return None

        all_keys = ["dfact_meV", "dfactprime_meV", "v0", "b0_GPa", "b1"]
        if what is None:
            keys = all_keys
        else:
            what = list_strings(what)
            if what[0].startswith("-"):
                # Exclude keys
                what = [w[1:] for w in what]
                keys = [k for k in all_keys if k not in what]
            else:
                keys = what

        # Get reference entry
        reference = df_database(xc=xc).get_entry(symbol=self.symbol, code=code)
        print("Reference data:", reference)

        # Get DataFrame.
        frame = self.get_pdframe(trial, *keys)
        ecuts = np.array(frame["ecut"])

        import matplotlib.pyplot as plt
        if ax_list is None:
            fig, ax_list = plt.subplots(nrows=len(keys), ncols=1, sharex=True, squeeze=False)
            ax_list = ax_list.ravel()
        else:
            fig = plt.gcf()

        if len(keys) != len(ax_list):
            raise ValueError("len(keys)=%s != len(ax_list)=%s" % (len(keys), len(ax_list)))

        for i, (ax, key) in enumerate(zip(ax_list, keys)):
            values = np.array(frame[key])
            refval = getattr(reference, key)
            # Plot difference pseudo - ref.
            #print("ecuts", ecuts, "values", values)
            psmae_diff = values - refval
            ax.plot(ecuts, psmae_diff, "o-")

            # Add vertical lines at hints.
            if self.has_hints:
                vmin, vmax = psmae_diff.min(), psmae_diff.max()
                for acc in self.ALL_ACCURACIES:
                    ax.vlines(self["hints"][acc]["ecut"], vmin, vmax,
                              colors=self.ACC2COLOR[acc], linestyles="dashed")

            ax.grid(True)
            ax.set_ylabel(r"$\Delta$" + key)
            if i == len(keys) - 1: ax.set_xlabel("Ecut [Ha]")

            xmin, xmax = min(ecuts), max(ecuts)
            if key == "dfactprime_meV":
                # Add horizontal lines (used to find hints for ecut).
                last = values[-1]
                for pad, acc in zip(self.ATOLS, self.ALL_ACCURACIES):
                    color = self.ACC2COLOR[acc]
                    ax.hlines(y=last + pad, xmin=xmin, xmax=xmax, colors=color, linewidth=1.5, linestyles='dashed')
                    ax.hlines(y=last - pad, xmin=xmin, xmax=xmax, colors=color, linewidth=1.5, linestyles='dashed')
                # Set proper limits so that we focus on the relevant region.
                #ax.set_ylim(last - 1.1*self.ATOLS[0], last + 1.1*self.ATOLS[0])
            else:
                ax.hlines(y=0., xmin=xmin, xmax=xmax, colors="black", linewidth=2, linestyles='dashed')

        plt.tight_layout()

        return fig

    @add_fig_kwargs
    def plot_gbrv_eos(self, struct_type, ax=None, with_soc=False, **kwargs):
        """
        Uses Matplotlib to plot the EOS computed with the GBRV setup

        Args:
            ax: matplotlib :class:`Axes` or None if a new figure should be created.
            with_soc: If True, the results obtained with SOC are plotted (if available).

        ================  ==============================================================
        kwargs            Meaning
        ================  ==============================================================
        cmap              Color map. default `jet`
        ================  ==============================================================

        Returns:
            `matplotlib` figure or None if the GBRV test is not present
        """
        trial = "gbrv_" + struct_type
        if with_soc: trial += "_soc"

        # Handle missing entries: noble gases, Hg ...
        if trial not in self:
            cprint("dojo report does not contain trial: %s" % str(trial), "red")
            return None

        ax, fig, plt = get_ax_fig_plt(ax)
        cmap = kwargs.pop("cmap", plt.get_cmap("jet"))

        # Get DataFrame.
        frame = self.get_pdframe(trial, "volumes", "etotals")
        ecuts = frame["ecut"]

        for i, ecut in enumerate(ecuts):
            # Subframe with this value of ecut.
            ecut_frame = frame.loc[frame["ecut"] == ecut]
            assert ecut_frame.shape[0] == 1
            # Extract volumes and energies for this ecut.
            volumes = (np.array(list(ecut_frame["volumes"].values), dtype=np.float)).flatten()
            etotals = (np.array(list(ecut_frame["etotals"].values), dtype=np.float)).flatten()

            eos_fit = EOS.Quadratic().fit(volumes, etotals)
            label = "ecut %.1f" % ecut if i % 2 == 0 else ""
            label = "ecut %.1f" % ecut
            eos_fit.plot(ax=ax, text=False, label=label, color=cmap(i/len(ecuts), alpha=1), show=False)

        plt.tight_layout()

        return fig

    @add_fig_kwargs
    def plot_gbrv_convergence(self, ax_list=None, with_soc=False, **kwargs):
        """
        Uses Matplotlib to plot the convergence of the GBRV parameters wrt ecut.

        Args:
            ax_list: List of matplotlib Axes, if ax_list is None a new figure is created
            with_soc: If True, the results obtained with SOC are plotted (if available).

        Returns:
            `matplotlib` figure. None if the GBRV test is not present.
        """
        import matplotlib.pyplot as plt
        stypes = ("fcc", "bcc")

        for stype in stypes:
            trial = "gbrv_" + stype
            if with_soc: trial += "_soc"
            if trial not in self:
                cprint("dojo report does not contain trial: %s" % str(trial), "red")
                return None

        if ax_list is None:
            fig, ax_list = plt.subplots(nrows=len(stypes), ncols=1, sharex=True, squeeze=False)
            ax_list = ax_list.ravel()
        else:
            fig = plt.gcf()

        if len(stypes) != len(ax_list):
            raise ValueError("len(stypes)=%s != len(ax_list)=%s" % (len(stypes), len(ax_list)))

        for i, (ax, stype) in enumerate(zip(ax_list, stypes)):
            trial = "gbrv_" + stype
            if with_soc: trial += "_soc"

            # Get DataFrame
            frame = self.get_pdframe(trial, "a0_rel_err")
            ecuts = np.array(frame["ecut"])
            values = np.array(frame["a0_rel_err"])

            ax.grid(True)
            ax.set_ylabel(r"$\Delta$" + trial + "a0_rel_err")

            # Plot difference pseudo - ref.
            ax.plot(ecuts, values, "bo-")
            ax.hlines(y=0., xmin=min(ecuts), xmax=max(ecuts), colors="black", linewidth=2, linestyles='dashed')
            if i == len(ax_list) - 1: ax.set_xlabel("Ecut [Ha]")

            # Add vertical lines at hints.
            if self.has_hints:
                vmin, vmax = values.min(), values.max()
                for acc in self.ALL_ACCURACIES:
                    ax.vlines(self["hints"][acc]["ecut"], vmin, vmax,
                              colors=self.ACC2COLOR[acc], linestyles="dashed")

        plt.tight_layout()

        return fig

    @add_fig_kwargs
    def plot_phonon_convergence(self, ax_list=None, with_soc=False, **kwargs):
        """
        Plot the convergence of the phonon modes wrt ecut.

        Args:
            ax_list: List of matplotlib Axes, if ax_list is None a new figure is created
            with_soc: If True, the results obtained with SOC are plotted (if available).

        Returns:
            `matplotlib` figure. None if the GBRV test is not present.
        """
        trial = "phgamma" if not with_soc else "phgamma_soc"
        if trial not in self:
            cprint("dojo report does not contain trial: %s" % str(trial), "red")
            return None

        frame = self.get_pdframe(trial, "asr2_phfreqs_mev", "noasr_phfreqs_mev")
        ecuts = np.array(frame["ecut"])
        num_modes = len(frame["asr2_phfreqs_mev"][0])

        # Build array with frequencies computed at the different ecut.
        asr2_phfreqs = np.empty((num_modes, len(ecuts)))
        noasr_phfreqs = np.empty((num_modes, len(ecuts)))

        mev2cmm1 = 8.065
        for ie, ecut in enumerate(ecuts):
            # Subframe with this value of ecut.
            ecut_frame = frame.loc[frame["ecut"] == ecut]
            asr2_phfreqs[:, ie] = mev2cmm1 * np.array(list(ecut_frame["asr2_phfreqs_mev"].values))
            noasr_phfreqs[:, ie] = mev2cmm1 * np.array(list(ecut_frame["noasr_phfreqs_mev"].values))

        import matplotlib.pyplot as plt
        from cycler import cycler
        fig, ax_list = plt.subplots(nrows=4, sharex=True)

        for ax in ax_list:
            ax.grid(True)
            ax.set_prop_cycle(cycler('color', ['r', 'g', 'b', 'y']) + cycler('linestyle', ['-', '--', ':', '-.']))

        for mu in range(num_modes):
            phecut = asr2_phfreqs[mu]
            ax_list[0].plot(ecuts, phecut, "o-")
            ax_list[0].set_ylabel(r"$\omega$ (asr=2)")
            values = phecut - phecut[-1]
            ax_list[1].plot(ecuts, values, "o-")
            ax_list[1].set_ylabel(r"$\omega-\omega_{max}$")

            # Add vertical lines at hints.
            if self.has_hints:
                vmin, vmax = values.min(), values.max()
                for acc in self.ALL_ACCURACIES:
                    ax_list[1].vlines(self["hints"][acc]["ecut"], vmin, vmax,
                                      colors=self.ACC2COLOR[acc], linestyles="dashed")

            phecut = noasr_phfreqs[mu]
            ax_list[2].plot(ecuts, phecut, "o-")
            ax_list[2].set_ylabel(r"$\omega$ (noasr)")
            ax_list[3].plot(ecuts, phecut - phecut[-1], "o-")
            ax_list[3].set_ylabel(r"$\omega-\omega_{max}$")

        # Adjust limits.
        fact = 0.05
        phmin, phmax = asr2_phfreqs.min(), asr2_phfreqs.max()
        if phmin == 0.0: phmin = -1
        ax_list[0].set_ylim(phmin - fact * abs(phmin), phmax + fact * abs(phmax))
        phmin, phmax = noasr_phfreqs.min(), noasr_phfreqs.max()
        if phmin == 0.0: phmin = -1
        ax_list[2].set_ylim(phmin - fact * abs(phmin), phmax + fact * abs(phmax))
        ax_list[-1].set_xlabel("Ecut [Ha]")

        fig.suptitle("Phonons in cm-1")

        plt.tight_layout()

        return fig

    @add_fig_kwargs
    def plot_ebands(self, ecut=None, with_soc=False, **kwargs):
        """
        Plot electronic band structure.

        Args:
            ecut:
            with_soc: If True, the results obtained with SOC are plotted (if available).

        ================  =============================
        kwargs            Meaning
        ================  =============================
        width             Gaussian broadening in eV
        step              Step of the DOS mesh in eV
        ================  =============================

        Returns:
            `matplotlib` figure. None if the ebands test is not present.
        """
        trial = "ghosts" if not with_soc else "ghosts_soc"
        if trial not in self:
            cprint("dojo report does not contain trial: %s" % str(trial), "red")
            return None
        ecut = list(self[trial].keys())[-1]

        from abipy.electrons.ebands import ElectronBands
        ebands = ElectronBands.from_dict(self[trial][ecut]["ebands"])
        edos = ebands.get_edos(width=kwargs.pop("width", 0.3), step=kwargs.pop("step", 0.1))

        # Try to detect possible ghost states by looking at the dispersion of the bands.
        dless_states = ebands.dispersionless_states(deltae=0.05, kfact=0.9)
        if not dless_states:
            print("No dispersionless state detected")
        else:
            print("Found %s dispersionless states" % len(dless_states))
            for i, s in enumerate(dless_states):
                print("[%d]" % i, s)

        return ebands.plot_with_edos(edos, show=False, **kwargs)

    @add_fig_kwargs
    def plot_raren_convergence(self, xc, with_soc=False, plot_diffs=False, **kwargs):
        """
        Plot the convergence of the total energy wrt ecut using (etotal obtained with the initial value
        of the lattice paramenter used to start the structural relaxation) as well as the convergence of the
        lattice parameter wrt to ecut.

        Args:
            xc: String or XcFunc object specifying the XC functional. E.g "PBE" or XcFunc.from_name("PBE")
            with_soc: If True, the results obtained with SOC are plotted (if available).
            plot_diffs: plot the difference with respect to the last value is stead of absolute values

        Returns:
            `matplotlib` figure. None if the ebands test is not present.
        """
        trial = "raren_relax" if not with_soc else "raren_relax_soc"
        if trial not in self:
            cprint("dojo report does not contain trial: %s" % str(trial), "red")
            return None

        # Energy is in eV/atom.
        key2ylabel = {"initial_energy_ev_per_atom": r"$\Delta E$ [meV/natom]", "relaxed_a": "$a$ [Angstrom]"}
        keys = list(key2ylabel.keys())
        data = self.get_pdframe(trial, *keys)

        import matplotlib.pyplot as plt
        fig, ax_list = plt.subplots(nrows=len(keys), ncols=1, sharex=True, squeeze=False)
        ax_list = ax_list.ravel()

        table = raren_database(xc).table

        ecuts = np.array(data["ecut"])
        xmin, xmax = min(ecuts), max(ecuts)
        for i, (ax, key) in enumerate(zip(ax_list, keys)):
            values = np.array(data[key])
            if key == "initial_energy_ev_per_atom":
                values = values * 1000
            diffs = values - values[-1]
            if plot_diffs:
                ax.plot(ecuts, diffs, "o-")
                vmin, vmax = diffs.min(), diffs.max()
                prec_list = [0.01, 0.005, 0.001] if key == "relaxed_a" else [1000, 200, 10]
                for prec in prec_list:
                    ax.hlines(y=prec, xmin=xmin, xmax=xmax, linewidth=1.5, linestyles='dashed')
                    ax.hlines(y=-prec, xmin=xmin, xmax=xmax, linewidth=1.5, linestyles='dashed')
            else:
                ax.plot(ecuts, values, "o-")
                vmin, vmax = values.min(), values.max()
                if key == "relaxed_a":
                    y = table["ref"][self.symbol]
                    ax.hlines(y=y, xmin=xmin, xmax=xmax, colors="b", linewidth=1.5, linestyles='dashed')
            # Add vertical lines at hints.
            if self.has_hints:
                for acc in self.ALL_ACCURACIES:
                    ax.vlines(self["hints"][acc]["ecut"], vmin, vmax, colors=self.ACC2COLOR[acc], linestyles="dashed")

            ax.grid(True)
            ax.set_ylabel(key2ylabel[key])
            if i == len(keys) - 1: ax.set_xlabel("Ecut [Ha]")

        plt.tight_layout()

        return fig

    #def get_raren_dataframe(self):
    #    db = raren_database(self.xc)
    #    return db.table.loc[self.symbol]


######################
## Pandas DataFrame ##
######################

class DojoDataFrame(pd.DataFrame):
    """
    Extends pandas `DataFrame` adding helper functions.
    """
    # For each trial, the quantities that will be stored in the dataframe.
    _TRIALS2KEY = {
        "ecut": "ecut",
        "deltafactor": "dfact_meV",
        "gbrv_bcc": "gbrv_bcc_a0_rel_err",
        "gbrv_fcc": "gbrv_fcc_a0_rel_err",
        "phgamma": "asr2_phfreqs_mev"
        #"phgamma": ["asr2_phfreqs_mev", "noasr_phfreqs_mev"]
    }

    _TRIALS2YLABEL = {
        "ecut": "Ecut [Ha]",
        "deltafactor": r"$\Delta$-factor [meV]",
        "gbrv_bcc": r"BCC $\Delta a_0$ (%)",
        "gbrv_fcc": r"FCC $\Delta a_0$ (%)",
        "phgamma": "asr2_phfreqs_mev"
    }

    # The frame has its own list so that one can easily change the
    # entries that should be analyzed by modifying this attributes.
    #ALL_TRIALS = DojoReport.ALL_TRIALS
    ALL_TRIALS = list(_TRIALS2KEY.keys())

    ALL_ACCURACIES = DojoReport.ALL_ACCURACIES

    ACC2PLTOPTS = dict(
        low=dict(color="orange"),
        normal=dict(color="green"),
        high=dict(color="red"),
    )

    #ACC2COLOR = {k: ACC2PLTOPTS[k]["color"] for k in ACC2PLTOPTS}

    for v in ACC2PLTOPTS.values():
        v.update(linewidth=2, linestyle='dashed', marker='o', markersize=8)
    del v

    @classmethod
    def from_json_file(cls, path, **kwargs):
        """
        Read the object from the json file `path`.
        kwargs are passed to `pandas.read_json`.
        """
        new = pd.read_json(path_or_buf=path, **kwargs)
        new.__class__ = cls
        return new

    @classmethod
    def from_pseudos(cls, pseudos):
        """
        Buid a pandas :class:`DataFrame` with the most important parameters
        extracted from the `DOJO_REPORT` section of each pseudo in the table.

        Returns: (frame, errors)

        where frame is the pandas :class:`DataFrame` and errors is a list of errors
        encountered while trying to read the `DOJO_REPORT` from the pseudopotential file.
        """
        accuracies = ["low", "normal", "high"]

        # For each trial, the quantities that will be stored in the dataframe.
        trial2keys = {
            "deltafactor": ["dfact_meV", "dfactprime_meV", "v0", "b0_GPa", "b1"],
            "gbrv_bcc": ["a0_rel_err"],
            "gbrv_fcc": ["a0_rel_err"],
            "phgamma": ["asr2_phfreqs_mev", "noasr_phfreqs_mev"]
        }

        rows, names, errors = [], [], []
        eapp = errors.append

        for p in pseudos:
            names.append(p.basename)
            d = {"symbol": p.symbol, "Z": p.Z, "Z_val": p.Z_val, "l_max": p.l_max,
                 "nlcc": p.has_nlcc, 'lmax': p.l_max} #"filepath": p.filepath}

            if not p.has_dojo_report:
                eapp("Cannot find dojo_report in %s" % p.basename)
                continue
            report = p.dojo_report
            d.update({"validated": report.isvalidated})

            ecut_acc = {}
            # read hints
            for acc in accuracies:
                try:
                    d.update({acc + "_ecut_hint": report['hints'][acc]['ecut']})
                    ecut_acc[acc] = report['hints'][acc]['ecut']
                except KeyError:
                    # using -1 for non existing values facilitates plotting
                    d.update({acc + "_ecut_hint": -1.0})
                    ecut_acc[acc] = -1

            for acc in accuracies:
                d[acc + "_ecut"] = ecut_acc[acc]

            try:
                for trial, keys in trial2keys.items():
                    data = report.get(trial, None)
                    if data is None:
                        eapp("No %s for %s" % (trial, p.basename))
                        continue

                    # if the current trial has an entry for this ecut change nothing, else we take the
                    # smallest, the middle and the highest ecut available for this trials
                    # precausion, normally either there are hints or not. in the second case they are all set to -1
                    ecut_acc_trial = dict(
                        low=sorted(data.keys())[0],
                        normal=sorted(data.keys())[int(len(data.keys())/2)],
                        high=sorted(data.keys())[-1],
                    )

                    for acc in accuracies:
                        d[acc + "_ecut"] = ecut_acc[acc]

                    for acc in accuracies:
                        ecut = ecut_acc[acc] if ecut_acc[acc] in data.keys() else ecut_acc_trial[acc]
                        # store the actual ecut for this trial
                        d.update({acc + "_ecut_" + trial: ecut})
                        if trial.startswith('ph'):
                            ecuts = data
                            d.update({acc + "_" + trial: data[ecut]})
                        else:
                            if trial.startswith("gbrv"):
                                d.update({acc + "_" + trial + "_" + k: float(data[ecut][k]) for k in keys})
                            else:
                                d.update({acc + "_" + k: float(data[ecut][k]) for k in keys})

            except Exception as exc:
                cprint("%s raised %s" % (p.basename, exc), "magenta")
                eapp("%s raised %s" % (p.basename, exc))

            rows.append(d)

        # Build sub-class of pandas.DataFrame
        return cls(rows, index=names), errors

    def myrows(self):
        """
        Return list with the row indices available in the dataframe.
        """
        rows = []
        for index, entry in self.iterrows():
            element = Element.from_Z(entry.Z)
            if element.row not in rows: rows.append(element.row)
        return sorted(rows)

    def myfamilies(self):
        """
        Return list of families available in the dataframe.
        """
        pd_families = [
           "noble_gas", "transition_metal", "rare_earth_metal", "metalloid",
           "alkali", "alkaline", "halogen", "lanthanoid", "actinoid",
        ]

        my_families = set()
        for index, entry in self.iterrows():
            element = Element.from_Z(entry.Z)
            for fam in pd_families:
                # e.g element.is_alkaline
                if getattr(element, "is_" + fam):
                    my_families.add(fam)
                    break

        return sorted(my_families)

    def select_rows(self, rows):
        """
        Select a list of rows of the periodic table.
        Rows can be a integer or list of integers.
        Return new :class:`DojoDataFrame`.
        """
        if not isinstance(rows, (list, tuple)): rows = [rows]
        rows = set(rows)

        data = []
        for index, entry in self.iterrows():
            element = Element.from_Z(entry.Z)
            if element.row in rows: data.append(entry)

        return self.__class__(data=data)

    def select_family(self, family):
        """
        Select a particular family of elements in the periodic table.
        Return new :class:`DojoDataFrame`.
        """
        data = []
        for index, entry in self.iterrows():
            element = Element.from_Z(entry.Z)
            # e.g element.is_alkaline
            if getattr(element, "is_" + family): data.append(entry)

        return self.__class__(data=data)

    def tabulate(self, columns=None, stream=sys.stdout):
        if columns is None:
            accuracies = self.ALL_ACCURACIES
            columns = [acc + "_dfact_meV" for acc in accuracies]
            columns += [acc + "_ecut" for acc in accuracies]
            columns += [acc + "_gbrv_fcc_a0_rel_err" for acc in accuracies]
            columns += [acc + "_gbrv_bcc_a0_rel_err" for acc in accuracies]

        stream.write(tabulate(self[columns], headers="keys", tablefmt="grid", floatfmt=".2f"))
        #return self[columns].to_html()

    def select_best(self):
        """
        Return dataframe with the best entries selected according to the deltafactor.
        """
        sortby, ascending = "high_dfact_meV", True
        rows, names = [], []
        for name, group in self.groupby("symbol"):
            # Sort group and select best pseudo depending on sortby and ascending.
            best = group.sort_values(sortby, ascending=ascending).iloc[0]
            names.append(name)
            #print(best.name, best.keys())
            l = {k: getattr(best, k) for k in ("name", "Z", "Z_val",
                                               "high_dfact_meV", "high_ecut_deltafactor",
                                               "high_gbrv_bcc_a0_rel_err", "high_gbrv_fcc_a0_rel_err"
                                              )}

            rows.append(l)

        best_frame = pd.DataFrame(rows, index=names)
        best_frame = best_frame.sort_values("Z")
        return best_frame

    ##################
    # Plotting tools #
    ##################

    @add_fig_kwargs
    def plot_hist(self, what="dfact_meV", bins=400, **kwargs):
        import matplotlib.pyplot as plt
        fig, ax_list = plt.subplots(nrows=len(self.ALL_ACCURACIES), ncols=1,
                                    sharex=True, sharey=False, squeeze=True)

        for acc, ax in zip(self.ALL_ACCURACIES, ax_list):
            col = acc + "_" + what
            self[col].plot(ax=ax, kind="bar", label=col)

        return fig

    @add_fig_kwargs
    def plot_trials(self, trials="all", accuracies="all", **kwargs):
        import matplotlib.pyplot as plt
        trials = self.ALL_TRIALS if trials == "all" else list_strings(trials)
        accuracies = self.ALL_ACCURACIES if accuracies == "all" else list_strings(accuracies)

        fig, ax_list = plt.subplots(nrows=len(trials), ncols=1, sharex=True, sharey=False, squeeze=True)

        # See also http://matplotlib.org/examples/pylab_examples/barchart_demo.html
        #print("frame keys:", self.keys())
        #print("phgamma:", self["normal_phgamma"])
        for i, (trial, ax) in enumerate(zip(trials, ax_list)):
            # FIXME
            if trial == "phgamma": continue
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
                start, end = 0.0, min(15, maxval)
                ax.set_ylim(start, end)
                #ax.yaxis.set_ticks(np.arange(start, end, 0.1))

            #if stepsize is not None:
            #    start, end = ax.get_ylim()
            #    ax.yaxis.set_ticks(np.arange(start, end, stepsize))

            plt.setp(ax.xaxis.get_majorticklabels(), rotation=25)

        return fig


class DfGbrvDataFrame(pd.DataFrame):
    """
    Extends pandas DataFrame adding helper functions.
    """
    @classmethod
    def from_pseudos(cls, pseudos, raise_if_none_dojoreport=False):
        """
        Build and return a pandas :class:`DataFrame` in the form.

            basename     deltafactor  df_prime  gbrv_bcc  gbrv_fcc  symbol   md5
            H-high.psp8  0.074830     1.258340  0.028904  0.024726  H        5863396c90149cbe12af496141bde0d0
            ...

        where `gbrv_bcc` and `gbrv_fcc` are the relative errors (in percentage) wrt the AE calculations.

        Args:
            raise_if_none_dojoreport: If True, a ValueError is raised if one of the pseudo does not
                have the dojo_report else a warning is emitted.
        """

        _TRIALS2KEY = {
            "deltafactor": "dfact_meV",
            "gbrv_bcc": "a0_rel_err",
            "gbrv_fcc": "a0_rel_err",
        }

        rows = []
        for p in pseudos:
            # Extract the dojo_report
            if not p.has_dojo_report:
                msg = "%s does not have the dojo_report" % p.filepath
                if not raise_if_none_dojoreport:
                    cprint(msg, "magenta")
                    continue
                else:
                    raise ValueError(msg)

            report = p.dojo_report
            row = dict(basename=p.basename, symbol=p.symbol, md5=p.md5)

            for trial, key in _TRIALS2KEY.items():
                # Get results as function of ecut
                try:
                    data = report[trial]
                except KeyError:
                    cprint("%s does not have %s" % (p.basename, trial), "red")
                    continue

                # Extract the value with highest ecut.
                high_ecut = list(data.keys())[-1]
                row.update({trial: data[high_ecut][key]})
                if trial == "deltafactor":
                    row.update(dict(df_prime=data[high_ecut]["dfactprime_meV"]))

            rows.append(row)

        return cls(rows)

    @add_fig_kwargs
    def plot_dfgbrv_dist(self, **kwargs):
        """
        Plot four distribution plots for the deltafactor, deltafactor prime and the
        relative errors for the GBRV fcc, bcc structures.

        Return: `matplotlib` figure.
        """
        import matplotlib.pyplot as plt
        fig, ax_list = plt.subplots(nrows=2, ncols=2, squeeze=True)
        ax_list = ax_list.ravel()

        import seaborn as sns
        for ax, col in zip(ax_list.ravel(), ["deltafactor", "gbrv_fcc", "df_prime", "gbrv_bcc"]):
            values = self[col].dropna()
            sns.distplot(values, ax=ax, rug=True, hist=True, kde=False, label=col, bins=kwargs.pop("bins", 50))

            # Add text with Mean or (MARE/RMSRE)
            text = []; app = text.append
            if col in ("deltafactor", "df_prime"):
                app("Mean = %.2f" % values.mean())
            else:
                app("MARE = %.2f" % values.abs().mean())
                app("RMSRE = %.2f" % np.sqrt((values**2).mean()))

            ax.text(0.8, 0.8, "\n".join(text), transform=ax.transAxes)

        return fig
