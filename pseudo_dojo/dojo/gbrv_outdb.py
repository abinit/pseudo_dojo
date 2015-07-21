# coding: utf-8
"""The GBRV results for binary and ternary compunds"""
from __future__ import division, print_function, unicode_literals

#import abc
import sys
import os
#import six
import json
import numpy as np

from collections import OrderedDict, MutableMapping, defaultdict
from warnings import warn
from monty.io import FileLock
from monty.string import list_strings
from atomicfile import AtomicFile
from pandas import DataFrame
from monty.collections import AttrDict, dict2namedtuple
from monty.functools import lazy_property
from pymatgen.core.periodic_table import sort_symbols_by_Z
#from pymatgen.core.units import Ha_to_eV
from pymatgen.util.plotting_utils import add_fig_kwargs, get_ax_fig_plt
from pseudo_dojo.core.pseudos import DojoTable
from pseudo_dojo.refdata.gbrv.database import gbrv_database, species_from_formula


import logging
logger = logging.getLogger(__name__)


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
        """Construct the object from a dictionary."""
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
                symbol = p.symbol
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
        # Validate input.
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
        Return a namedtuple with the paramenters to be used to run
        the job for the specified accuracy.
        """
        assert accuracy in self.ACCURACIES
        # TODO
        #@ecut = max(p.hint_for_accuracy(accuracy).ecut for p in self.pseudos)        
        ecut, pawecutdg = 6, None

        return dict2namedtuple(formula=self.formula, struct_type=self.struct_type, accuracy=accuracy, 
                               pseudos=self.pseudos, ecut=ecut, pawecutdg=pawecutdg)

    def compute_err(self, reference="ae", accuracy="normal"):
        """
        Return namedtuple with absolute and relative error.
        None if data is not available.
        """
        # Get the reference results
        gbrv_db = gbrv_database()
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
        plot the EOS computed with the deltafactor setup.

        Args:
            ax: matplotlib :class:`Axes` or None if a new figure should be created.
            
        Returns:
            `matplotlib` figure.
        """
        ax, fig, plt = get_ax_fig_plt(ax)

        if not self.has_data(accuracy): return fig
        d = self["accuracy"]

        num_sites, volumes, etotals = d["num_sites"], np.array(d["volumes"]), np.array(d["etotals"])
        from pymatgen.io.abinitio.eos import EOS
        eos = EOS.Quadratic()

        # Use same fit as the one employed for the deltafactor.
        eos_fit = eos.fit(volumes/num_sites, etotals/num_sites)

        label = "ecut %.1f" % d["ecut"]
        eos_fit.plot(ax=ax, text=False, label=label,show=False) # color=cmap(i/num_ecuts, alpha=1), 

        return fig


class GbrvOutdb(MutableMapping):
    """
    This object stores the results for the GBRV tests (binary and ternary compounds.
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
                new[formula] = [GbrvRecord.from_dict(d, new.struct_type, dojo_pptable) for d in dict_list]

            return new

    def __init__(self, dojo_dir, dojo_pptable):
        gbrv_db = gbrv_database()
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
        d = {}
        d["struct_type"] = self.struct_type
        d["dojo_dir"] = self.dojo_dir
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
        gbrv_db = gbrv_database()
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
        Reset all the records whose status is in status_list so that we can resubmit them.
        Return number of records that have been resetted.
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
        Build a pandas :class:`DataFrame` with the most important results.

        Args:
            reference:
            pptable: :class:`PseudoTable` object. If given. the frame will contains only the 
                entries with pseudopotential in pptable.

        Returns:
            frame: pandas :class:`DataFrame` 
        """
        rows = []

        for formula, records in self.items():
            for rec in records:
                d = {"formula": formula, "basenames": set(p.basename for p in rec.pseudos)}
                d.update({"md5": {p.symbol: p.md5 for p in rec.pseudos}})

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

    #ALL_ACCURACIES = ("normal", "high")

    @lazy_property
    def symbols(self):
        """List with the element symbols present in the table sorted by Z."""
        symbols = set()
        for idx, row in self.iterrows(): 
            symbols.update(row.md5.keys())

        return sort_symbols_by_Z(symbols)

    #@lazy_property
    #def pseudo_metas(self):
    #    d = {}
    #    for idx, row in self.iterrows(): 
    #        d.update(row.pseudo_metas)
    #    return d

    def pprint(self, **kwargs):
        """Pretty-print"""
        frame = self[["formula", "high_rel_err", "basenames"]]
        s = frame.to_string(index=False)
        print(s)
        print("")
        #print(frame.describe())
        #print("")

    def subframe_for_pseudo(self, pseudo, best_for_accuracy=None):
        """
        Extract the rows with the given pseudo. Return new `GbrvDataFrame`.

        Args:
            pseudo: :class:`Pseudo` object or string with the pseudo basename.
            best_for_accuracy: If not None, the returned frame will contain one
                entry for formula. This entry has the `best` relative error
                i.e. it's the one with the minimum absolute error.
        """
        pname = pseudo.basename if hasattr(pseudo, "basename") else pseudo
                                                                             
        # Extract the rows containing this pseudo.
        rows = [row for idx, row in self.iterrows() if pname in row.basenames]

        new = self.__class__(rows)
        if best_for_accuracy is None: 
            return new

        # Handle best_for_accuracy
        key = best_for_accuracy + "_rel_err"

        #groups = new.groupby("formula")
        ##groups = new.groupby("basenames").groups
        #for group in groups:
        #    print("group")
        #    print(group)
        #    #print(group[key])
        #raise NotImplementedError()

    def subframe_for_symbol(self, symbol):
        """Extract the rows with the given element symbol. Return new `GbrvDataFrame`."""
        rows = [row for idx, row in self.iterrows() if symbol in row.md5.keys()]
        return self.__class__(rows)

    @add_fig_kwargs
    def plot_error_pseudo(self, pseudo, ax=None, **kwargs):
        frame = self.subframe_for_pseudo(pseudo)

        ynames = ["normal_rel_err", "high_rel_err"]
        style = ["-o"] * len(ynames)

        ax, fig, plt = get_ax_fig_plt(ax)
        frame.plot("formula", ynames, style=style, grid=True, ax=ax
                        #kind="scatter"
                       )

        #xticks = ax.get_xticks()
        #ax.set_xticklabels(frame.formula)
        return fig

    @add_fig_kwargs
    def plot_allpseudos_with_symbol(self, symbol, accuracy="normal", **kwargs):
        # For each pseudo:
        # Extract the sub-frame for this pseudo and keep the rows with the 
        # best result for the given accuracy
        ax, fig, plt = get_ax_fig_plt(ax)
        key = accuracy + "_rel_err"

        # Find all pseudos with the given symbol in the table.
        frame = self.subframe_for_symbol(symbol)

        #for pseudo in frame.pseudos:
        #    frame.plot_error_pseudo(pseudo, ax=None)

        return fig

