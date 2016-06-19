# coding: utf-8
"""The GBRV results for binary and ternary compunds"""
from __future__ import division, print_function, unicode_literals

import os
import json
import numpy as np

from collections import OrderedDict, MutableMapping, defaultdict
from warnings import warn
from monty.io import FileLock
from monty.string import list_strings
from atomicfile import AtomicFile
from pandas import DataFrame
from monty.collections import dict2namedtuple #AttrDict,
from monty.functools import lazy_property
from pymatgen.core.periodic_table import Element #sort_symbols_by_Z
from pymatgen.util.plotting_utils import add_fig_kwargs, get_ax_fig_plt
from pymatgen.analysis.eos import EOS
from pseudo_dojo.core.pseudos import DojoTable
from pseudo_dojo.refdata.gbrv.database import gbrv_database, species_from_formula


import logging
logger = logging.getLogger(__name__)


def sort_symbols_by_Z(symbols):
    """
    Given a list of element symbols, sort the strings according to Z,
    Return sorted list.

    >>> assert sort_symbols_by_Z(["Si", "H"]) == ["H", "Si"]
    """
    return list(sorted(symbols, key=lambda s: Element(s).Z))


class GbrvRecord(dict):
    """
    Example of entry of LiCl:

        record = {
            formula: "LiCl",
            pseudos_metadata: {
                "Li": {basename: "Li-s-high.psp8", md5: 312},
                "Cl": {basename: "Cl.psp8", md5: 562}],
            }
            normal: {ecut:, a:, v0: , b0:, b1:},
            high: {ecut:, a:, v0: , b0:, b1},
        }

        where results is a dictionary

        "normal": {
            "ecut": 6,
            "v0": 31.72020565768123,
            "a0": 5.024952898489712,
            "etotals": [-593.3490598305451, ...],
            "b0": 4.148379951739942,
            "b1": Infinity,
            "volumes": [31.410526411353832, ...],
            "num_sites": 2
        }
    """
    ACCURACIES = ("normal", "high")

    STATUS_LIST = (None, "scheduled" ,"failed")

    @classmethod
    def from_dict(cls, d, struct_type, dojo_pptable):
        """
        Construct the object from a dictionary.
        """
        d = d.copy()
        new = cls(struct_type, d.pop("formula"), d.pop("pseudos_metadata"), dojo_pptable)

        for acc in cls.ACCURACIES:
            new[acc] = d.pop(acc)

        assert not d
        return new

    def as_dict(self):
        """Dict representation."""
        return {k: self[k] for k in self}

    def __init__(self, struct_type, formula, pseudos_or_dict, dojo_pptable):
        """
        Initialize the record for the chemical formula and the list of
        pseudopotentials.
        """
        keys = ("basename", "Z_val", "l_max", "md5")

        #if isinstance(pseudos_or_dict, (list, tuple)):
        if all(hasattr(p, "as_dict") for p in pseudos_or_dict):
            def get_info(p):
                """Extract the most important info from the pseudo."""
                #symbol = p.symbol
                d = p.as_dict()
                return {k: d[k] for k in keys}

            meta = {p.symbol: get_info(p) for p in pseudos_or_dict}
            pseudos = pseudos_or_dict

        else:
            meta = pseudos_or_dict
            for v in meta.values():
                assert set(v.keys()) == set(keys)

            def pmatch(ps, esymb, d):
                return (ps.md5 == d["md5"] and
                        ps.symbol == esymb and
                        ps.Z_val == d["Z_val"] and
                        ps.l_max == d["l_max"])

            pseudos = []
            for esymb, d in meta.items():
                for p in dojo_pptable.pseudo_with_symbol(esymb, allow_multi=True):
                    if pmatch(p, esymb, d):
                        pseudos.append(p)
                        break
                else:
                    raise ValueError("Cannot find pseudo:\n %s\n in dojo_pptable" % str(d))

        super(GbrvRecord, self).__init__(formula=formula, pseudos_metadata=meta,
                                         normal=None, high=None)

        self.pseudos = DojoTable.as_table(pseudos)
        self.dojo_pptable = dojo_pptable
        self.struct_type = struct_type

    def __eq__(self, other):
        """Two records can be compared for equality."""
        if other is None: return False
        if self["formula"] != other["formula"]: return False
        if len(self.pseudos) != len(other.pseudos): return False
        if self.struct_type != other.struct_type: return False

        for p in self.pseudos:
            try:
                o = other.pseudos.pseudo_with_symbol(p.symbol)
            except ValueError:
                return False

            if p != o: return False

        return True

    def __ne__(self, other):
        return not self == other

    def __str__(self):
        return json.dumps(self.as_dict(), indent=4, sort_keys=False)

    @property
    def formula(self):
        return self["formula"]

    def add_results(self, accuracy, results):
        # TODO Validate input.
        assert accuracy in self.ACCURACIES
        #assert set(data.keys()) == set(["ecut", "a0", "v0" , "b0", "b1"])

        if self[accuracy] is not None:
            logger.warning("Overwriting results for %s" % self.formula)

        self[accuracy] = results

    def has_data(self, accuracy):
        """True if the record contains computed data for this accuracy."""
        return self[accuracy] not in self.STATUS_LIST

    def matches_pseudos(self, pseudos):
        """Return True if the list of `pseudos` matches the one in the record."""
        d1 = {p.symbol: p for p in self.pseudos}
        d2 = {p.symbol: p for p in pseudos}
        return d1 == d2

    def get_jobparams(self, accuracy):
        """
        Return a namedtuple with the parameters to be used
        to run the job for the specified accuracy.
        """
        assert accuracy in self.ACCURACIES
        # TODO
        #@ecut = max(p.hint_for_accuracy(accuracy).ecut for p in self.pseudos)
        ecut, pawecutdg = 6, None

        return dict2namedtuple(formula=self.formula, struct_type=self.struct_type, accuracy=accuracy,
                               pseudos=self.pseudos, ecut=ecut, pawecutdg=pawecutdg)

    def compute_err(self, reference="ae", accuracy="normal"):
        """
        Return namedtuple with the absolute and the relative error.
        None if data is not available.
        """
        # Get the reference results
        gbrv_db = gbrv_database(xc=None) # FIXME
        gbrv_entry = gbrv_db.get_entry(self.formula, self.struct_type)
        ref_a = getattr(gbrv_entry, reference)

        # Get our value, Return None if not computed.
        if not self.has_data(accuracy): return None
        d = self[accuracy]

        our_a = d["a0"]
        abs_err = our_a - ref_a
        rel_err = 100 * abs_err / ref_a

        return dict2namedtuple(a0=our_a, abs_err=abs_err, rel_err=rel_err)

    @add_fig_kwargs
    def plot_eos(self, ax=None, accuracy="normal", **kwargs):
        """
        Plot the equation of state.

        Args:
            ax: matplotlib :class:`Axes` or None if a new figure should be created.

        Returns:
            `matplotlib` figure.
        """
        ax, fig, plt = get_ax_fig_plt(ax)

        if not self.has_data(accuracy): return fig
        d = self["accuracy"]

        num_sites, volumes, etotals = d["num_sites"], np.array(d["volumes"]), np.array(d["etotals"])

        # Perform quadratic fit.
        eos = EOS.Quadratic()
        eos_fit = eos.fit(volumes/num_sites, etotals/num_sites)

        label = "ecut %.1f" % d["ecut"]
        eos_fit.plot(ax=ax, text=False, label=label,show=False) # color=cmap(i/num_ecuts, alpha=1),

        return fig


class GbrvOutdb(MutableMapping):
    """
    This object stores the results for the GBRV tests (binary and ternary compounds).
    This is a base class.
    """
    struct_type = None

    @classmethod
    def new_from_dojodir(cls, dojo_dir):
        """
        Initialize the object from a top level directory that
        contains pseudopotentials in the PseudoDojo format.
        """
        # Construct the full table of pseudos from dojodir
        dojo_pptable = DojoTable.from_dojodir(dojo_dir)

        # Here I initialize the object with default data (None).
        new = cls(dojo_dir, dojo_pptable)

        for formula, species in new.gbrv_formula_and_species:

            # Find all the possible combinations of dojo_pptable compatible with this formula!
            comb_list = dojo_pptable.all_combinations_for_elements(set(species))
            #print("Found ",len(comb_list), " combinations for %s" % formula)

            # Add record for this formula.
            for pplist in comb_list:
                new[formula].append(GbrvRecord(new.struct_type, formula, pplist, dojo_pptable))

        return new

    @classmethod
    def from_file(cls, filepath):
        """
        Initalize the object from a file in json format.
        """
        with open(filepath, "rt") as fh:
            d = json.loads(fh.read())
            struct_type = d.pop("struct_type")
            dojo_dir = d.pop("dojo_dir")

            if cls.struct_type != struct_type:
                for subcls in cls.__subclasses__():
                    if subcls.struct_type == struct_type:
                        cls = subcls
                        break
                else:
                    raise ValueError("Cannot find subclass associated to %s" % struct_type)

            # Construct the full table of pseudos from dojodir
            dojo_pptable = DojoTable.from_dojodir(dojo_dir)

            # Here I initialize the object with the data read from file.
            new = cls(dojo_dir, dojo_pptable)
            for formula, dict_list in d.items():
                new[formula] = [GbrvRecord.from_dict(d, new.struct_type, dojo_pptable)
                                for d in dict_list]

            return new

    def __init__(self, dojo_dir, dojo_pptable):
        gbrv_db = gbrv_database(xc=None) # FIXME
        ord_keys = gbrv_db.tables[self.struct_type].keys()

        self.dojo_dir = os.path.basename(dojo_dir)
        self.dojo_pptable = dojo_pptable

        self.data = OrderedDict()
        for key in ord_keys:
            self[key] = []

    def __str__(self):
        #return str(self.data)
        return self.to_json()

    @property
    def filepath(self):
        # TODO: dojo_dir is relative path that should be converted
        dirpath = self.dojo_dir
        return os.path.join(dirpath, self.basename)

    def to_json(self):
        d = dict(struct_type=self.struct_type, dojo_dir=self.dojo_dir)
        for formula, records in self.items():
            d[formula] = [rec.as_dict() for rec in records]

        return json.dumps(d, indent=4, sort_keys=False)

    def json_write(self, filepath=None):
        """Write new file with locking mechanism."""
        # FIXME: This cannot be called my multiple processes!
        filepath = self.filepath if filepath is None else filepath
        with FileLock(filepath):
            with AtomicFile(filepath, mode="wt") as fh:
                fh.write(self.to_json())

    @classmethod
    def update_record(cls, filepath, formula, accuracy, pseudos, results):
        """Update database."""
        with FileLock(filepath):
            outdb = cls.from_file(filepath)

            rec = outdb.find_record(formula, pseudos)
            rec.add_results(accuracy, results)

            with AtomicFile(filepath, mode="wt") as fh:
                fh.write(outdb.to_json())

    # ABC protocol.
    def __iter__(self):
        return self.data.__iter__()

    def __len__(self):
        return self.data.__len__()

    def __getitem__(self, key):
        return self.data[key]

    def __setitem__(self, key, value):
        self.data[key] = value

    def __delitem__(self, key):
        raise ValueError("Cannot delete key %s" % key)
    # End ABC

    def find_record(self, formula, pseudos):
        """
        Find the record associated to the specified formula
        and the list of pseudos. Return None if record is not present.
        """
        if formula not in self: return None

        records = self[formula]
        try:
            i = records.index(GbrvRecord(self.struct_type, formula, pseudos, self.dojo_pptable))
        except ValueError:
            return None

        return records[i]

    def has_record(self, record):
        return record in self[record.formula]

    def find_jobs_torun(self, max_njobs=3, select_formulas=None):
        """
        Find entries whose results have not been yet calculated.

        Args:
            select_formulas:
        """
        jobs, got = [], 0
        #if max_njobs == -1: max_njobs = np.inf

        for formula, records in self.items():
            if got >= max_njobs: break
            if select_formulas is not None and formula not in select_formulas:
                continue

            for rec in records:
                for accuracy in rec.ACCURACIES:
                    data = rec[accuracy]
                    # TODO: Better treatment of failed!
                    if data in ("scheduled", "failed"): continue
                    if got < max_njobs and data is None:
                        got += 1
                        jobs.append(rec.get_jobparams(accuracy))
                        rec[accuracy] = "scheduled"

        # Update the database.
        if jobs:
            self.json_write()

        return jobs

    @lazy_property
    def gbrv_formula_and_species(self):
        """
        Return a list of tuples. Each tuple contains the formula and the list
        of chemical species in the same order as in formula.
        Example: ("LiF2", ["Li", "F", "F"])
        """
        gbrv_db = gbrv_database(xc=None) # FIXME
        formulas = gbrv_db.tables[self.struct_type].keys()

        items = []
        for formula in formulas:
            items.append((formula, species_from_formula(formula)))

        return items

    def check_update(self):
        """
        Check consistency between the pseudo potential table and the database and upgrade it
        This usually happens when new pseudopotentials have been added to the dojo directory.
        (very likely) or when pseudos have been removed (unlikely!)

        Returns: namedtuple with the following attributes.
            nrec_removed
            nrec_added
        """
        nrec_removed, nrec_added = 0, 0
        missing = defaultdict(list)

        for formula, species in self.gbrv_formula_and_species:
            # Get **all** the possible combinations for these species.
            comb_list = self.dojo_pptable.all_combinations_for_elements(set(species))

            # Check consistency between records and pseudos!
            # This is gonna be slow if we have several possibilities!
            records = self[formula]
            recidx_found = []
            for pseudos in comb_list:
                for i, rec in enumerate(records):
                    if rec.matches_pseudos(pseudos):
                        recidx_found.append(i)
                        break

                else:
                    missing[formula].append(pseudos)

            # Remove stale records (if any)
            num_found = len(recidx_found)
            if  num_found != len(records):
                num_stale = len(records) - num_found
                print("Found %s stale records" % num_stale)
                nrec_removed += num_stale
                self[formula] = [records[i] for i in recidx_found]

        if missing:
            for formula, pplist in missing.items():
                for pseudos in pplist:
                    nrec_removed += 1
                    self[formula].append(GbrvRecord(self.struct_type, formula, pseudos, self.dojo_pptable))

        if missing or nrec_removed:
            print("Updating database.")
            self.json_write()

        return dict2namedtuple(nrec_removed=nrec_removed, nrec_added=nrec_added)

    def reset(self, status_list="failed"):
        """
        Reset all the records whose status is in `status_list` so that we can resubmit them.
        Return the number of records that have been resetted.
        """
        status_list = list_strings(status_list)
        count = 0

        for formula, records in self.items():
            for rec in records:
                for accuracy in rec.ACCURACIES:
                    if rec[accuracy] in status_list:
                        count += 1
                        rec[accuracy] = None

        return count

    @add_fig_kwargs
    def plot_errors(self, reference="ae", accuracy="normal", ax=None, **kwargs):
        """
        Plot the error wrt the reference values.

        Args:
            ax: matplotlib :class:`Axes` or None if a new figure should be created.

        Returns:
            `matplotlib` figure
        """
        ax, fig, plt = get_ax_fig_plt(ax)

        ax.grid(True)
        ax.set_xlabel('r [Bohr]')

        xs, ys_abs, ys_rel = [], [], []

        for formula, records in self.items():
            for rec in records:
                e = rec.compute_err(reference=reference, accuracy=accuracy)
                if e is None: continue
                xs.append(rec.formula)
                ys_abs.append(e.abs_err)
                ys_rel.append(e.rel_err)

        if not xs:
            warn("No entry available for plotting")

        ax.scatter(range(len(ys_rel)), ys_rel, s=20) #, c='b', marker='o', cmap=None, norm=None)
        #ax.scatter(xs, ys_rel, s=20) #, c='b', marker='o', cmap=None, norm=None)

        return fig

    def get_dataframe(self, reference="ae", pptable=None, **kwargs):
        """
        Build and return a :class:`GbrvDataFrame` with the most important results.

        Args:
            reference:
            pptable: :class:`PseudoTable` object. If given. the frame will contains only the
                entries with pseudopotential in pptable.

        Returns:
            frame: pandas :class:`DataFrame`
        """
        def get_df(p):
            dfact_meV, df_prime = None, None

            if p.has_dojo_report:
                try:
                    data = p.dojo_report["deltafactor"]
                    high_ecut = list(data.keys())[-1]
                    dfact_meV = data[high_ecut]["dfact_meV"]
                    df_prime = data[high_ecut]["dfactprime_meV"]
                except KeyError:
                    pass

            return dict(dfact_meV=dfact_meV, df_prime=df_prime)

        def get_meta(p):
            """Return dict with pseudo metadata."""
            meta = {"basename": p.basename, "md5": p.md5}
            meta.update(get_df(p))
            return meta

        rows = []
        for formula, records in self.items():
            for rec in records:
                d = dict(formula=formula, struct_type=self.struct_type,
                         basenames=set(p.basename for p in rec.pseudos),
                         pseudos_meta={p.symbol: get_meta(p) for p in rec.pseudos},
                         symbols={p.symbol for p in rec.pseudos},
                        )

                has_data = 0
                for acc in ("normal", "high"):
                    e = rec.compute_err(reference=reference, accuracy=acc)
                    if e is None: continue
                    has_data += 1
                    d.update({acc + "_a0": e.a0, acc + "_rel_err": e.rel_err, #acc + "_abs_err": e.abs_err,
                             })

                if has_data:
                    rows.append(d)

        # Build sub-class of pandas.DataFrame
        new = GbrvDataFrame(rows)
        if pptable is None: return new

        raise NotImplementedError()
        # Remove rows that are not in pptable.
        #rows = []
        #for index, entry in self.iterrows():
        #    md5 = entry.md5
        #    if not pptable.has_md5_signatures(md5): continue
        #    rows.append(entry)

        #return GbrvDataFrame(rows)


class RocksaltOutdb(GbrvOutdb):
    """Results for the rocksalt structures."""
    struct_type = "rocksalt"
    basename = struct_type + ".json"


class PeroviskiteOutdb(GbrvOutdb):
    """Results for the ABO3 structures."""
    struct_type = "ABO3"
    basename = struct_type + ".json"


class HalfHeuslersOutdb(GbrvOutdb):
    """Results for the half-Heuslers structures."""
    struct_type = "hH"
    basename = struct_type + ".json"


class GbrvDataFrame(DataFrame):
    """
    Extends pandas :class:`DataFrame` by adding helper functions.

    The frame has the structure:

          a0         high_rel_err  normal_rel_err     basenames
    CsCl  7.074637           NaN       -0.188528      set(Cs_basename, Cl_basename)
    BeO   3.584316           NaN       -1.799555      {...}
    """
    ALL_ACCURACIES = ("normal", "high")

    @classmethod
    def from_dojodir(cls, dojo_dir, exclude_basenames=None):
        """
        Initialize the object from a top level directory that
        contains pseudopotentials in the PseudoDojo format.

        Args:
            exclude_basenames: Optional string or list of strings with the
                pseudo basenames to be excluded.
        """
        # Construct the full table of pseudos from dojodir
        dojo_pptable = DojoTable.from_dojodir(dojo_dir, exclude_basenames=exclude_basenames)
        return cls.from_dojotable(dojo_pptable)

    @classmethod
    def from_dojotable(cls, table):

        _TRIALS2KEY = {
            #"deltafactor": "dfact_meV",
            "gbrv_bcc": "a0_rel_err",
            "gbrv_fcc": "a0_rel_err",
        }

        rows = []

        # This is for scatter.
        hack = False
        if hack:
            xs, ys, sizes = [], [], []

        for p in table:
            # Extract GBRV results from dojo_report
            report = p.dojo_report

            # Add column with deltafactor prime.
            # FIXME should depend on the accuracy.
            try:
                df_dict = report["deltafactor"]
            except KeyError:
                print("%s does not have deltafactor" % p.basename)
                continue

            df_ecuts, df_values = df_dict.keys(), df_dict.values()
            last_ecut = list(sorted(df_ecuts))[-1]
            df_prime = df_dict[last_ecut]["dfactprime_meV"]
            df = df_dict[last_ecut]["dfact_meV"]
            print(p.basename, "last_ecut:", last_ecut, "df", df, "df_prime", df_prime)

            for trial in _TRIALS2KEY:
                # Get results as function of ecut
                try:
                    data = report[trial]
                except KeyError:
                    print("%s does not have %s" % (p.basename, trial))
                    continue

                struct_type = {"gbrv_bcc": "bcc", "gbrv_fcc": "fcc"}[trial]
                ecuts = data.keys()
                try:
                    values = np.array([float(data[ecut]["a0_rel_err"]) for ecut in ecuts])
                except KeyError:
                    print("Problem in a0_rel_err with %s" % p.basename)
                    continue

                row = dict(formula=trial, struct_type=struct_type,
                           basenames=set([p.basename]),
                           pseudos_meta={p.symbol: {"basename": p.basename, "md5": p.md5}},
                           symbols={p.symbol})

                for acc in cls.ALL_ACCURACIES:
                    # FIXME: ecut should depend on accuracy.
                    # for the time being we get the last one
                    ecut = ecuts[-1]

                    row.update({acc + "_rel_err": data[ecut]["a0_rel_err"]})
                    row.update({acc + "_df": df, acc + "_df_prime": df_prime})

                rows.append(row)

                if hack:
                    # Hack for scatter
                    if trial == "gbrv_fcc":
                        if "gbrv_bcc" in report:
                            xs.append(data[ecut]["a0_rel_err"])
                    else:
                        if "gbrv_fcc" in report:
                            ys.append(data[ecut]["a0_rel_err"])

                    if "gbrv_bcc" in report and "gbrv_fcc" in report:
                        #sizes.append(df)
                        sizes.append(df_prime)

        # Scatter hack
        if hack:
            ax, fig, plt = get_ax_fig_plt()
            sizes = np.array(sizes)
            sizes = 60 * sizes / sizes.max()
            sizes = sizes**2

            ax.scatter(xs, ys, s=sizes, alpha=0.5) #c=close,
            ax.grid(True)

            l = np.linspace(np.min(xs), np.max(xs), num=50)
            ax.plot(l, l)

            #ax.set_xlabel("bcc rel error")
            #ax.set_ylabel("fcc_rel_error")
            fig.tight_layout()
            plt.show()

        new = cls(rows)
        return new

    @lazy_property
    def symbols(self):
        """List with the element symbols present in the table sorted by Z."""
        symbols = set()
        for idx, row in self.iterrows():
            symbols.update(row.symbols)

        return sort_symbols_by_Z(symbols)

    #def multiple_pseudos(self):
    #    # Loop over the rows. Collect all the {symbol: md5}
    #    meta_set, multiple = set(), []
    #    for idx, row in self.iterrows():
    #        for esymb, meta in row.pseudos_meta.items():
    #            meta = tuple([(k, v) for k, v in meta.items()])
    #            if meta not in meta_set:
    #                meta_set.add(meta)
    #            else:
    #                multiple.append(meta)
    #    return [{k: v for k, v in t} for t in multiple]

    def print_info(self, **kwargs):
        """Pretty-print"""
        frame = self[["formula", "normal_rel_err", "high_rel_err", "basenames"]] # "symbols",
        s = frame.to_string(index=False)
        print(s)
        print("")
        #print(frame.describe())
        #print("")

        for col in ["normal_rel_err", "high_rel_err"]:
            print("For column:" , col)
            print("mean(abs)", self[col].abs().mean())
            print("RMS:", np.sqrt((self[col]**2).sum() / len(self)))

    def select_badguys(self, accuracy="high", atol=0.5):
        col = accuracy + "_rel_err"
        #bad = [row for idx, row in self.iterrows() if abs(row[col]) > atol]
        #new = self.__class__(bad)
        new = self[self[col].abs() > atol]
        new.__class__ = self.__class__
        return new

    @add_fig_kwargs
    def plot_hist(self, ax=None, **kwargs):
        """Histogram plot."""
        ax, fig, plt = get_ax_fig_plt(ax)
        import seaborn as sns

        ax.grid(True)
        #for acc in ("normal", "high"):
        for acc in ("high",):
            col = acc + "_rel_err"
            values = self[col].dropna()
            sns.distplot(values, ax=ax, rug=True, hist=False, label=col)

            # Add text with Mean or (MARE/RMSRE)
            text = []; app = text.append
            app("MARE = %.2f" % values.abs().mean())
            app("RMSRE = %.2f" % np.sqrt((values**2).mean()))

            ax.text(0.8, 0.8, "\n".join(text), transform=ax.transAxes)

        return fig

    def subframe_for_symbol(self, symbol):
        """Extract the rows with the given element symbol. Return new `GbrvDataFrame`."""
        # Extract the rows containing this pseudo and create new frame.
        rows = []
        for idx, row in self.iterrows():
            if symbol not in row.symbols: continue
            meta = row.pseudos_meta[symbol]

            pseudo_basename = meta["basename"]
            dfact_meV, df_prime = meta["dfact_meV"], meta["df_prime"]
            row = row.set_value("dfact_meV", dfact_meV)
            row = row.set_value("dfactprime_meV", df_prime)
            row = row.set_value("pseudo_basename", pseudo_basename)

            rows.append(row)

        return self.__class__(rows)

    def subframe_for_pseudo(self, pseudo, best_for_acc=None):
        """
        Extract the rows with the given pseudo. Return new `GbrvDataFrame`.

        Args:
            pseudo: :class:`Pseudo` object or string with the pseudo basename.
            best_for_acc: If not None, the returned frame will contain one
                entry for formula. This entry has the `best` relative error
                i.e. it's the one with the minimum absolute error.
        """
        pname = pseudo.basename if hasattr(pseudo, "basename") else pseudo

        rows = []
        for idx, row in self.iterrows():
            if pname not in row.basenames: continue
            meta = row.pseudos_meta[p.symbol]
            row.set_value("pseudo_basename", pname)
            # Add values of deltafactor
            #dfact_meV, df_prime = extract_df(row)
            #row.set_value("dfact_meV", dfact_meV)
            #row.set_value("dfactprime_meV", df_prime)
            rows.append(row)

        new = self.__class__(rows)
        if best_for_acc is None:
            return new

        # Handle best_for_acc
        key = best_for_acc + "_rel_err"

        #raise NotImplementedError()
        #rows = []
        #for formula, group in new.groupby("formula"):
        #    iloc = group[key].abs().idxmin()
        #    #best = group.iloc(loc)
        #    best = new.iloc(iloc)
        #    print(type(best), best)
        #    #print(group[best])
        #    #rows.append(best.set_value("formula", formula))

        d = defaultdict(list)
        for idx, row in new.iterrows():
            d[row.formula].append((row, row[key]))

        rows = []
        for formula, values in d.items():
            best_row = sorted(values, key= lambda t: abs(t[1])) [0][0]
            rows.append(best_row)

        return self.__class__(rows)

    @add_fig_kwargs
    def plot_error_pseudo(self, pseudo, ax=None, **kwargs):
        #frame = self.subframe_for_pseudo(pseudo)
        frame = self.subframe_for_pseudo(pseudo, best_for_acc="high")

        ax, fig, plt = get_ax_fig_plt(ax)
        for acc in self.ALL_ACCURACIES:
            yname = acc + "_rel_err"
            frame.plot("formula", yname, grid=True, ax=ax, style="-o", label=acc)

        ax.set_ylim(-1.0, +1.0)
        ax.legend(loc="best")
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=70)

    @add_fig_kwargs
    def boxplot(self, ax=None, **kwargs):
        import seaborn as sns
        ax, fig, plt = get_ax_fig_plt(ax)

        ax = sns.boxplot(self["high_rel_err"], groupby=self.formula) #, orient="h")

        plt.setp(ax.xaxis.get_majorticklabels(), rotation=70)
        ax.set_ylim(-0.5, 0.5)
        ax.grid(True)

        return fig

    @add_fig_kwargs
    def stripplot_symbol(self, symbol, **kwargs):
        frame = self.subframe_for_symbol(symbol)

        import seaborn as sns
        sns.set(style="whitegrid", palette="pastel")
        #ax, fig, plt = get_ax_fig_plt(None)

        import matplotlib.pyplot as plt
        fig, ax_list = plt.subplots(nrows=2, ncols=1, squeeze=True)
        ax0, ax1 = ax_list.ravel()

        ax1 = ax_list[1]
        sns.stripplot(x="pseudo_basename", y="high_rel_err", data=frame, hue="formula", ax=ax1,
                      jitter=True, size=10, marker="o", edgecolor="gray", alpha=.25, #palette="Set2",
        )

        ax1.grid(True)
        ax1.axhline(y=-0.4, linewidth=2, color='r', linestyle="--")
        ax1.axhline(y=0.0, linewidth=2, color='k', linestyle="--")
        ax1.axhline(y=+0.4, linewidth=2, color='r', linestyle="--")

        # Plot the deltafactor for the different pseudos on another Axes.
        xlabels = ax1.xaxis.get_majorticklabels()

        #xs, ys, ls  = [], [], []
        rows = []
        for xlabel in xlabels:
            #print(dir(xlabel))
            (x, y), basename = xlabel.get_position(), xlabel.get_text()

            g = frame[frame["pseudo_basename"] == basename]
            df = g.iloc[0]["dfact_meV"]
            #print(x, y, df)
            rows.append(dict(pseudo_basename=basename, dfact_meV=df))
            #xs.append(x)
            #ys.append(df)
            #ls.append(basename)

        #print(xs)
        #ax0.plot(xs, ys, "-o")
        frame = DataFrame(rows)
        frame.plot("pseudo_basename", "dfact_meV", ax=ax0, style="-o")

        return fig

    #@add_fig_kwargs
    #def plot_allpseudos_with_symbol(self, symbol, accuracy="normal", **kwargs):
    #    # For each pseudo:
    #    # Extract the sub-frame for this pseudo and keep the rows with the
    #    # best result for the given accuracy
    #    ax, fig, plt = get_ax_fig_plt(None)
    #    key = accuracy + "_rel_err"
    #    # Find all pseudos with the given symbol in the table.
    #    frame = self.subframe_for_symbol(symbol)
    #    #for pseudo in frame.pseudos:
    #    #    frame.plot_error_pseudo(pseudo, ax=None)
    #    #import seaborn as sns
    #    # Initialize a grid of plots with an Axes for each walk
    #    #grid = sns.FacetGrid(df, col="walk", hue="walk", col_wrap=5, size=1.5)
    #    #grid = sns.FacetGrid(frame, col=key)#, hue="walk", col_wrap=5, size=1.5)
    #    # Draw a horizontal line to show the starting point
    #    #grid.map(plt.axhline, y=0, ls=":", c=".5")
    #    # Draw a line plot to show the trajectory of each random walk
    #    #grid.map(plt.plot, "formula", key, marker="o", ms=4)
    #    # Adjust the tick positions and labels
    #    #grid.set(xticks=np.arange(5), yticks=[-3, 3],
    #    #         xlim=(-.5, 4.5), ylim=(-3.5, 3.5))
    #    # Adjust the arrangement of the plots
    #    #grid.fig.tight_layout(w_pad=1)
    #    return fig

    @add_fig_kwargs
    def hist_allpseudos_with_symbols(self, symbol, ax=None, **kwargs):
        import seaborn as sns
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        # Find all entries with this symbol and add new column with the basename
        frame = self.subframe_for_symbol(symbol)
        #frame["pseudo_name"] = [entry.pseudos_meta[symbol]["basename"] for index, entry in frame.iterrows()]

        # Group by basename and plot.
        grouped = frame.groupby("pseudo_basename")

        for name, group in grouped:
            print(name) #; print(group)
            acc = "high"
            col = acc + "_rel_err"
            s = group[col].dropna()
            if len(s) in [0, 1]: continue
            print(s)
            sns.distplot(s, ax=ax, rug=True, hist=True, kde=False, label=name)
            ax.set_xlim(-0.5, 0.5)

        ax.axvline(x=0, linewidth=2, color='k', linestyle="--")
        ax.axvline(x=0.2, linewidth=2, color='r', linestyle="--")
        ax.axvline(x=-0.2, linewidth=2, color='r', linestyle="--")
        ax.legend(loc="best")

        return fig

    @add_fig_kwargs
    def plot_formulas(self, ax=None, **kwargs):
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        ynames = ["normal_rel_err", "high_rel_err"]
        ynames = ["high_rel_err"]

        self.plot("formula", ynames, ax=ax, style="o", grid=True)

        plt.setp(ax.xaxis.get_majorticklabels(), rotation=70)
        ax.set_ylim(-1, 1)

        return fig

    @add_fig_kwargs
    def join_plot(self, **kwargs):
        import seaborn as sns
        sns.set(style="darkgrid", color_codes=True)

        ax, fig, plt = get_ax_fig_plt()

        #tips = sns.load_dataset("tips")
        #g = sns.jointplot("total_bill", "tip", data=tips, kind="reg",
        #                  xlim=(0, 60), ylim=(0, 12), color="r", size=7)
        #print_full_frame(self[["formula", "high_df_prime", "basenames"]])

        newcol = "abs(high_rel_err)"

        self[newcol] = self["high_rel_err"].abs()

        g = sns.jointplot("high_df_prime", "abs(high_rel_err)", data=self, kind="reg",)
                          #xlim=(0, 60), ylim=(0, 12), color="r", size=7)

        g = sns.jointplot("high_df", "abs(high_rel_err)", data=self, kind="reg",)
                          #xlim=(0, 60), ylim=(0, 12), color="r", size=7)

        # Remove the column
        self.drop([newcol], axis=1)
        return fig

    #def scatter(self, **kwarags):
    #    ax, fig, plt = get_ax_fig_plt()
    #    xs
    #    ys
    #    #ax.scatter(delta1[:-1], delta1[1:], c=close, s=volume, alpha=0.5)
    #    return fig
