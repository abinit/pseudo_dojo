# coding: utf-8
"""The GBRV results for binary and ternary compunds"""
from __future__ import division, print_function, unicode_literals

import os
import json
import numpy as np

from collections import OrderedDict, defaultdict
from monty.io import FileLock
from monty.string import list_strings
from atomicfile import AtomicFile
from pandas import DataFrame
from monty.collections import dict2namedtuple
from monty.functools import lazy_property
from pymatgen.core.periodic_table import Element
from pymatgen.util.plotting_utils import add_fig_kwargs, get_ax_fig_plt
from pymatgen.analysis.eos import EOS
from pseudo_dojo.core.pseudos import DojoTable, OfficialDojoTable
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

    "rocksalt": {
        "LiCl" = {
            formula: "LiCl",
            pseudos_metadata: {
                "Li": {basename: "Li-s-high.psp8", md5: 312},
                "Cl": {basename: "Cl.psp8", md5: 562}],
            }
            normal: results,
            high: results,
        },
    }

    where results is the dictionary:

        {"ecut": 6,
         "v0": 31.72020565768123,
         "a0": 5.024952898489712,
         "b0": 4.148379951739942,
         "num_sites": 2
         "etotals": [-593.3490598305451, ...],
         "volumes": [31.410526411353832, ...],
        }
    """

    ACCURACIES = ("normal", "high")
    STATUS_LIST = (None, "scheduled" ,"failed")

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

    #def __str__(self):
    #    return json.dumps(self.as_dict(), indent=4, sort_keys=False)

    #@property
    #def formula(self):
    #    return self["formula"]

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
        eos_fit.plot(ax=ax, text=False, label=label, show=False) # color=cmap(i/num_ecuts, alpha=1),
        return fig


class GbrvOutdb(dict):
    """
    Stores the results for the GBRV tests (binary and ternary compounds).
    This object is usually created via the class methods:

        GbrvOutdb.from_file and GbrvOutdb.new_from_table.
    """
    # The structures stored in the database.
    struct_types = ["rocksalt", "ABO3", "hH"]

    # The name of the json database should start with prefix.
    prefix = "gbrv_compouds_"

    @classmethod
    def new_from_table(cls, table, djson_path):
        """
        Initialize the object from an :class:`OfficialDojoTable` and a djson file.
        """
        djson_path = os.path.abspath(djson_path)
        dirname = os.path.dirname(djson_path)
        new = cls(path=os.path.join(dirname, cls.prefix + os.path.basename(djson_path)),
                  djson_path=djson_path, xc_name=table.xc.name)

        # Init subdictionaries e.g. {'ABO3': {'KHgF3': None, 'SbNCa3': None}, "hH": {...}}
        # These dictionaries will be filled afterwards with the GBRV results.
        gbrv = gbrv_database(xc=table.xc)
        for struct_type in cls.struct_types:
            new[struct_type] = {k: None for k in gbrv.tables[struct_type]}

        # Get a reference to the table.
        new.table = table
        return new

    @classmethod
    def from_file(cls, filepath):
        """
        Initalize the object from a file in json format.
        """
        with open(filepath, "rt") as fh:
            d = json.load(fh)
            new = cls(**d)
            #new["xc"] = new["xc"]
            #print("keys", new.keys())

        # Construct the full table of pseudos
        djpath = new["djson_path"]
        from pseudo_dojo.pseudos import dojotable_absdir
        if not os.path.exists(djpath):
            head, base = os.path.split(djpath)
            _, dirname = os.path.split(head)
            dojodir = dojotable_absdir(dirname)
            djpath = os.path.join(dojodir, base)

        new.table = OfficialDojoTable.from_djson_file(djpath)
        return new

    def iter_struct_formula_data(self):
        """Iterate over (struct_type, formula, data)."""
        for struct_type in self.struct_types:
            for formula, data in self[struct_type].items():
                yield struct_type, formula, data

    @property
    def path(self):
        return self["path"]

    def json_write(self, filepath=None):
        """
        Write data to file in JSON format.
        If filepath is None, self.path is used.
        """
        filepath = self.path if filepath is None else filepath
        with open(filepath, "wt") as fh:
            json.dump(self, fh, indent=-1, sort_keys=True)

    @classmethod
    def insert_results(cls, filepath, struct_type, formula, accuracy, pseudos, results):
        """
        Update the entry in the database.
        """
        with FileLock(filepath):
            outdb = cls.from_file(filepath)
            old_dict = outdb[struct_type][formula]
            if not isinstance(old_dict, dict): old_dict = {}
            old_dict[accuracy] = results
            outdb[struct_type][formula] = old_dict
            with AtomicFile(filepath, mode="wt") as fh:
                json.dump(outdb, fh, indent=-1, sort_keys=True) #, cls=MontyEncoder)

    def find_job_torun(self):
        """
        Find entries whose results have not been yet calculated.

        Args:
            select_formulas:
        """
        job = None
        for struct_type, formula, data in self.iter_struct_formula_data():
            if data in ("scheduled", "failed"): continue
            if data is None:
                symbols = list(set(species_from_formula(formula)))
                pseudos = self.table.pseudos_with_symbols(symbols)

                job = dict2namedtuple(formula=formula, struct_type=struct_type, pseudos=pseudos)
                self[struct_type][formula] = "scheduled"
                break

        # Update the database.
        if job is not None: self.json_write()
        return job

    def get_record(self, struct_type, formula):
        """
        Find the record associated to the specified structure type and chemical formula.
        Return None if record is not present.
        """
        d = self.get(struct_type)
        if d is None: return None
        data = d.get(formula)
        if data is None: return None

        raise NotImplementedError()
        #return GbrvRecord.from_data(data, struct_type, formula, pseudos)

    # TODO
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

    def reset(self, status_list="failed", write=True):
        """
        Reset all the records whose status is in `status_list` so that we can resubmit them.
        Return the number of records that have been resetted.
        """
        status_list = list_strings(status_list)
        count = 0
        for struct_type, formula, data in self.iter_struct_formula_data():
            if data in status_list:
                self[struct_type][formula] = None
                count += 1

        # Update the database.
        if count: self.json_write()
        return count

    #############################
    ### Post-processing tools ###
    #############################

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
        #ax.set_xlabel('r [Bohr]')

        xs, ys_abs, ys_rel = [], [], []

        for formula, records in self.items():
            for rec in records:
                e = rec.compute_err(reference=reference, accuracy=accuracy)
                if e is None: continue
                xs.append(rec.formula)
                ys_abs.append(e.abs_err)
                ys_rel.append(e.rel_err)

        if not xs:
            print("No entry available for plotting")
            return None

        ax.scatter(range(len(ys_rel)), ys_rel, s=20) #, c='b', marker='o', cmap=None, norm=None)
        #ax.scatter(xs, ys_rel, s=20) #, c='b', marker='o', cmap=None, norm=None)
        return fig

    def get_pdframe(self, reference="ae", pptable=None, **kwargs):
        """
        Build and return a :class:`GbrvCompoundDataFrame` with the most important results.

        Args:
            reference:
            pptable: :class:`PseudoTable` object. If given. the frame will contains only the
                entries with pseudopotential in pptable.

        Returns:
            frame: pandas :class:`DataFrame`
        """
        #def get_df(p):
        #    dfact_meV, df_prime = None, None
        #    if p.has_dojo_report:
        #        try:
        #            data = p.dojo_report["deltafactor"]
        #            high_ecut = list(data.keys())[-1]
        #            dfact_meV = data[high_ecut]["dfact_meV"]
        #            df_prime = data[high_ecut]["dfactprime_meV"]
        #        except KeyError:
        #            pass
        #    return dict(dfact_meV=dfact_meV, df_prime=df_prime)

        #def get_meta(p):
        #    """Return dict with pseudo metadata."""
        #    meta = {"basename": p.basename, "md5": p.md5}
        #    meta.update(get_df(p))
        #    return meta

        gbrv = gbrv_database(xc=self["xc_name"])

        rows, miss = [], 0
        for struct_type, formula, data in self.iter_struct_formula_data():
            if not isinstance(data, dict):
                miss += 1
                continue

            a0 = data["high"]["a0"]

            row = dict(formula=formula, struct_type=struct_type, this=a0,
                     #basenames=set(p.basename for p in rec.pseudos),
                     #pseudos_meta={p.symbol: get_meta(p) for p in rec.pseudos},
                     #symbols={p.symbol for p in rec.pseudos},
                     )

            # Add results from the database.
            entry = gbrv.get_entry(formula, struct_type)
            #row.update({code: getattr(entry, code) for code in entry.code_names})
            row.update({code: getattr(entry, code) for code in ["ae", "gbrv_paw"]})
            rows.append(row)

        if miss:
            print("There are %d missing entries in %s\n" % (miss, self.path))

        # Build sub-class of pandas.DataFrame and add relative error wrt AE results.
        frame = GbrvCompoundDataFrame(rows)
        frame["rel_err"] = 100 * (frame["this"] - frame["ae"]) / frame["ae"]
        return frame


class GbrvCompoundDataFrame(DataFrame):
    """
    Extends pandas :class:`DataFrame` by adding helper functions.

    The frame has the structure:

              ae formula  gbrv_paw  gbrv_uspp  pslib struct_type      this   vasp
        0   4.073     AlN     4.079      4.079  4.068    rocksalt  4.068777  4.070
        1   5.161    LiCl     5.153      5.151  5.160    rocksalt  5.150304  5.150
        2   4.911      YN     4.909      4.908  4.915    rocksalt  4.906262  4.906

    TODO: column with pseudos?

    basenames
    set(Cs_basename, Cl_basename)
    {...}
    """
    ALL_ACCURACIES = ("normal", "high")

    #def select_badguys(self, rtol=0.4):
        #new = self[abs(100 * (self["this"] - self["ae"]) / self["ae"]) > rtol].copy()
        #new["rel_err"] = 100 * (self["this"] - self["ae"]) / self["ae"]
        #new.__class__ = self.__class__
        #return new

    #def print_info(self, **kwargs):
    #    """Pretty-print"""
    #    frame = self[["formula", "normal_rel_err", "high_rel_err", "basenames"]] # "symbols",
    #    s = frame.to_string(index=False)
    #    print(s)
    #    print("")
    #    #print(frame.describe())
    #    #print("")

    #    for col in ["normal_rel_err", "high_rel_err"]:
    #        print("For column:" , col)
    #        print("mean(abs)", self[col].abs().mean())
    #        print("RMS:", np.sqrt((self[col]**2).sum() / len(self)))

    #def subframe_for_symbol(self, symbol):
    #    """Extract the rows with the given element symbol. Return new `GbrvCompoundDataFrame`."""
    #    # Extract the rows containing this pseudo and create new frame.
    #    rows = []
    #    for idx, row in self.iterrows():
    #        if symbol not in row.symbols: continue
    #        meta = row.pseudos_meta[symbol]

    #        pseudo_basename = meta["basename"]
    #        dfact_meV, df_prime = meta["dfact_meV"], meta["df_prime"]
    #        row = row.set_value("dfact_meV", dfact_meV)
    #        row = row.set_value("dfactprime_meV", df_prime)
    #        row = row.set_value("pseudo_basename", pseudo_basename)

    #        rows.append(row)

    #    return self.__class__(rows)

    @add_fig_kwargs
    def plot_errors_for_structure(self, struct_type, ax=None, **kwargs):
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        data = self[self["struct_type"] == struct_type].copy()

        colnames = ["this", "gbrv_paw"]
        for col in colnames:
            data[col + "_rel_err"] = 100 * (data[col] - data["ae"]) / data["ae"]
            data.plot(x="formula", y=col + "_rel_err", ax=ax, style="o-", grid=True)

        #print(data)
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

    def subframe_for_pseudo(self, pseudo, best_for_acc=None):
        """
        Extract the rows with the given pseudo. Return new `GbrvCompoundDataFrame`.

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
        if best_for_acc is None: return new

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
