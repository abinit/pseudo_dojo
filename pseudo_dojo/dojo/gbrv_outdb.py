# coding: utf-8
"""The GBRV results for binary and ternary compunds"""
from __future__ import division, print_function, unicode_literals

import os
import json
import numpy as np

from collections import OrderedDict, defaultdict, Counter
from monty.io import FileLock
from monty.string import list_strings
from atomicfile import AtomicFile
from pandas import DataFrame
from monty.collections import dict2namedtuple
from monty.functools import lazy_property
from monty.os.path import which
from pymatgen.core.periodic_table import Element
from pymatgen.util.plotting_utils import add_fig_kwargs, get_ax_fig_plt
from pymatgen.analysis.eos import EOS
from pseudo_dojo.core.pseudos import DojoTable, OfficialDojoTable
from pseudo_dojo.refdata.gbrv.database import gbrv_database, gbrv_code_names, species_from_formula
from pseudo_dojo.pseudos import as_dojo_path


import logging
logger = logging.getLogger(__name__)


def sort_symbols_by_Z(symbols):
    """
    Given a list of element symbols, sort the strings according to Z,
    Return sorted list.

    >>> assert sort_symbols_by_Z(["Si", "H"]) == ["H", "Si"]
    """
    return list(sorted(symbols, key=lambda s: Element(s).Z))


def print_frame(x):
    import pandas as pd
    with pd.option_context('display.max_rows', len(x),
                           'display.max_columns', len(list(x.keys()))):
        print(x)


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
        # Translate djson_path into path insides pseudos
        djpath = as_dojo_path(new["djson_path"])

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

    def find_jobs_torun(self, max_njobs):
        """
        Find entries whose results have not been yet calculated.

        Args:
            select_formulas:
        """
        jobs, got = [], 0
        for struct_type, formula, data in self.iter_struct_formula_data():
            if got == max_njobs: break
            if data in ("scheduled", "failed"): continue
            if data is None:
                symbols = list(set(species_from_formula(formula)))
                pseudos = self.table.pseudos_with_symbols(symbols)

                job = dict2namedtuple(formula=formula, struct_type=struct_type, pseudos=pseudos)
                self[struct_type][formula] = "scheduled"
                jobs.append(job)
                got += 1

        # Update the database.
        if jobs: self.json_write()
        return jobs

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
        if count and write: self.json_write()
        return count

    #############################
    ### Post-processing tools ###
    #############################

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
            row.update({code: getattr(entry, code) for code in gbrv_code_names})
            #row.update({code: getattr(entry, code) for code in ["ae", "gbrv_paw"]})
            rows.append(row)

        if miss:
            print("There are %d missing entries in %s\n" % (miss, self.path))

        # Build sub-class of pandas.DataFrame and add relative error wrt AE results.
        frame = GbrvCompoundDataFrame(rows)
        frame["rel_err"] = 100 * (frame["this"] - frame["ae"]) / frame["ae"]
        return frame

    def make_open_notebook(self, nbpath=None):
        """
        Generate an ipython notebook open it in the browser. Return system exit code.

        Raise:
            RuntimeError if jupyther or ipython are not in $PATH
        """
        path = self.write_notebook(nbpath=nbpath)

        if which("jupyter") is not None:
            return os.system("jupyter notebook %s" % path)

        if which("ipython") is not None:
            return os.system("ipython notebook %s" % path)

        raise RuntimeError("Cannot find neither jupyther nor ipython. Install them with `pip install`")

    def write_notebook(self, nbpath):
        """
        Write an ipython notebook. If `nbpath` is None, a temporay file is created.

        Returns:
            The path to the ipython notebook.

    See also:
        http://nbviewer.jupyter.org/github/maxalbert/auto-exec-notebook/blob/master/how-to-programmatically-generate-and-execute-an-ipython-notebook.ipynb
        """
        frame = self.get_pdframe()

        import nbformat
        nbv = nbformat.v4

        nb = nbv.new_notebook()
        nb.cells.extend([
            nbv.new_markdown_cell("# This is an auto-generated notebook"),
            nbv.new_code_cell("""\
from __future__ import print_function, division, unicode_literals
from IPython.display import display
import seaborn
%matplotlib notebook"""),

            nbv.new_code_cell("""\
from pseudo_dojo.dojo.gbrv_outdb import GbrvOutdb
outdb = GbrvOutdb.from_file('%s')
frame = outdb.get_pdframe()
frame.print_summary()""" % as_dojo_path(self.path)),

            nbv.new_code_cell("display(frame)"),
        ])

        for struct_type in frame.struct_types():
            nb.cells += [
                nbv.new_markdown_cell("## GBRV results for structure %s:" % struct_type),
                nbv.new_code_cell("fig = frame.plot_errors_for_structure('%s')" % struct_type),
                nbv.new_code_cell("fig = frame.plot_hist('%s')" % struct_type),
        ]

        nb.cells += [
            nbv.new_markdown_cell("## GBRV Compounds: relative errors as function of chemical element"),
            nbv.new_code_cell("fig = frame.plot_errors_for_elements()"),
            nbv.new_markdown_cell("## Bad guys:"),
            nbv.new_code_cell("bad, count = frame.select_bad_guys(reltol=0.4)"),
            nbv.new_code_cell("display(bad)"),
            nbv.new_code_cell("print(count)"),
        ]

        import io, tempfile
        if nbpath is None:
            _, nbpath = tempfile.mkstemp(suffix='.ipynb', text=True)

        with io.open(nbpath, 'wt', encoding="utf8") as f:
            nbformat.write(nb, f)

        return nbpath


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

    def struct_types(self):
        """List of structure types available in the dataframe."""
        return sorted(set(self["struct_type"]))

    def code_names(self):
        """List of code names available in the dataframe."""
        codes = sorted(c for c in self.keys() if c in gbrv_code_names)
        # Add this
        return ["this"] + codes

    def print_summary(self, choice="rms", **kwargs):
        """Print table with the rms of the different codes."""
        index, rows = [], []
        ref_code = "ae"
        for struct_type in self.struct_types():
            new = self[self["struct_type"] == struct_type].copy()
            index.append(struct_type)
            row = {}
            for code in self.code_names():
                if code == ref_code: continue
                values = (100 * (new[code] - new[ref_code]) / new[ref_code]).dropna()
                if choice == "rms":
                    row[code] = np.sqrt((values**2).mean())
                elif choice == "mare":
                    row[code] = values.abs().mean()
                else:
                    raise ValueError("Wrong value of choice: %s" % choice)
            rows.append(row)

        frame = DataFrame(rows, index=index)
        print(frame)

    def select_bad_guys(self, reltol=0.4):
        new = self[abs(100 * (self["this"] - self["ae"]) / self["ae"]) > reltol].copy()
        new["rel_err"] = 100 * (self["this"] - self["ae"]) / self["ae"]
        new.__class__ = self.__class__
        count = Counter()
        for idx, row in new.iterrows():
            for symbol in set(species_from_formula(row.formula)):
                count[symbol] += 1

        return new, count

    def remove_bad_guys(self, reltol=0.4):
        new = self[abs(100 * (self["this"] - self["ae"]) / self["ae"]) <= reltol].copy()
        new.__class__ = self.__class__
        new["rel_err"] = 100 * (self["this"] - self["ae"]) / self["ae"]
        return new

    def select_symbol(self, symbol):
        """
        Extract the rows whose formula contain the given element symbol.
        Return new `GbrvCompoundDataFrame`.
        """
        rows = []
        for idx, row in self.iterrows():
            if symbol not in species_from_formula(row.formula): continue
            #meta = row.pseudos_meta[symbol]
            #pseudo_basename = meta["basename"]
            #dfact_meV, df_prime = meta["dfact_meV"], meta["df_prime"]
            #row = row.set_value("dfact_meV", dfact_meV)
            #row = row.set_value("dfactprime_meV", df_prime)
            #row = row.set_value("pseudo_basename", pseudo_basename)
            rows.append(row)

        return self.__class__(rows)

    @add_fig_kwargs
    def plot_errors_for_structure(self, struct_type, ax=None, **kwargs):
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        data = self[self["struct_type"] == struct_type].copy()
        if not len(data):
            print("No results available for struct_type:", struct_type)
            return None

        colnames = ["this", "gbrv_paw"]
        for col in colnames:
            data[col + "_rel_err"] = 100 * (data[col] - data["ae"]) / data["ae"]
            #data[col + "_rel_err"] = abs(100 * (data[col] - data["ae"]) / data["ae"])
            data.plot(x="formula", y=col + "_rel_err", ax=ax, style="o-", grid=True)

        ax.set_title(struct_type)
        ax.set_ylabel("relative error %")
        # Set xticks and lables (show'em all)
        #print(data)
        ax.set_xticks(list(range(len(data.index))))
        ax.set_xticklabels(data["formula"], rotation="vertical")
        ax.tick_params(which='both', direction='out')
        ax.set_ylim(-1, 1)

        return fig

    @add_fig_kwargs
    def plot_hist(self, struct_type, ax=None, **kwargs):
        """
        Histogram plot.
        """
        #if codes is None: codes = ["ae"]
        ax, fig, plt = get_ax_fig_plt(ax)
        import seaborn as sns

        codes = ["this", "gbrv_paw"] #, "gbrv_uspp", "pslib", "vasp"]
        new = self[self["struct_type"] == struct_type].copy()
        for code in codes:
            values = (100 * (new[code] - new["ae"]) / new["ae"]).dropna()
            sns.distplot(values, ax=ax, rug=True, hist=False, label=code)

            if code == "this":
                # Add text with Mean or (MARE/RMSRE)
                text = []; app = text.append
                app("MARE = %.2f" % values.abs().mean())
                app("RMSRE = %.2f" % np.sqrt((values**2).mean()))
                ax.text(0.8, 0.8, "\n".join(text), transform=ax.transAxes)

        ax.grid(True)
        ax.set_xlabel("relative error %")
        ax.set_xlim(-0.8, 0.8)

        return fig

    @add_fig_kwargs
    def plot_errors_for_elements(self, ax=None, **kwargs):
        """
        Plot the relative errors associated to the chemical elements.
        """
        dict_list = []
        for idx, row in self.iterrows():
            rerr = 100 * (row["this"] - row["ae"]) / row["ae"]
            for symbol in set(species_from_formula(row.formula)):
                dict_list.append(dict(
                    element=symbol,
                    rerr=rerr,
                    formula=row.formula,
                    struct_type=row.struct_type,
                    ))

        frame = DataFrame(dict_list)
        order = sort_symbols_by_Z(set(frame["element"]))
        #print_frame(frame)

        import seaborn as sns
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        # Draw violinplot
        #sns.violinplot(x="element", y="rerr", order=order, data=frame, ax=ax, orient="v")

        # Box plot
        ax = sns.boxplot(x="element", y="rerr", data=frame, ax=ax, order=order, whis=np.inf, color="c")
        # Add in points to show each observation
        sns.stripplot(x="element", y="rerr", data=frame, ax=ax, order=order,
                      jitter=True, size=5, color=".3", linewidth=0)

        sns.despine(left=True)
        ax.set_ylabel("Relative error %")
        ax.grid(True)
        return fig
