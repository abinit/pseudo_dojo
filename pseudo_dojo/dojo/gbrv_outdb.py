# coding: utf-8
"""The GBRV results for binary and ternary compunds"""
from __future__ import division, print_function, unicode_literals

#import abc
#import sys
import os
#import six
import json
import numpy as np

from collections import OrderedDict, MutableMapping, defaultdict
from warnings import warn
from monty.io import FileLock
from monty.collections import AttrDict, dict2namedtuple
from monty.functools import lazy_property
#from monty.pprint import pprint_table
#from pymatgen.core.units import Ha_to_eV
from pymatgen.util.plotting_utils import add_fig_kwargs, get_ax_fig_plt
#from pymatgen.io.abinitio.pseudos import Pseudo
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
    """
    ACCURACIES = ("normal", "high")

    @classmethod
    def from_dict(cls, d, dojo_pptable):
        d = d.copy()
        new = cls(d.pop("formula"), d.pop("pseudos_metadata"), dojo_pptable)

        for acc in cls.ACCURACIES:
            new[acc] = d.pop(acc)

        assert not d
        return new

    def as_dict(self):
        return {k: self[k] for k in self}

    def __init__(self, formula, pseudos_or_dict, dojo_pptable):
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

    def __eq__(self, other):
        if other is None: return False
        if self["formula"] != other["formula"]: return False
        if len(self.pseudos) != len(other.pseudos): return False

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
        assert set(data.keys()) == set(["ecut", "a", "v0" , "b0", "b1"])

        self[accuracy] = results

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
        ecut = 6
        pawecutdg = None
        struct_type = None #?

        return dict2namedtuple(formula=self.formula, struct_type=struct_type, accuracy=accuracy, 
                               pseudos=self.pseudos, ecut=ecut, pawecutdg=pawecutdg)

    def compute_err(self, reference="ae", accuracy="normal"):
        # Get the reference results
        gbrv_db = gbrv_database()
        gbrv_entry = gbrv_db.get_entry(self.symbol)
        ref_a = getattr(gbrv_entry, reference)

        # Get our value, Return None if not computed.
        try:
            our_a = self[accuracy]["a"]
        except KeyError:
            return None

        abs_err = our_a - ref_a
        rel_err = 100 * abs_err / ref_a
        return dict2namedtuple(abs_err=abs_err, rel_err=rel_err)


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
                new[formula].append(GbrvRecord(formula, pplist, dojo_pptable))

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
                new[formula] = [GbrvRecord.from_dict(d, dojo_pptable) for d in dict_list]

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
    def basename(self):
        return self.struct_type + ".json"

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

        #d = self.data.copy()
        return json.dumps(d, indent=4, sort_keys=False)

    def json_write(self, filepath=None):
        """Write new file with locking mechanism."""
        filepath = self.filepath if filepath is None else filepath
        with FileLock(filepath):
            with open(filepath, "wt") as fh:
                fh.write(self.to_json())

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
            i = records.index(GbrvRecord(formula, pseudos, self.dojo_pptable))
        except ValueError:
            return None

        return records[i]

    def has_record(self, record):
        return record in self[record.formula]

    def find_jobs_torun(self, max_njobs=3):
        """
        Find entries whose results have not been yet calculated.
        """
        jobs, got = [], 0
        #if max_njobs == -1: max_njobs = np.inf

        for formula, records in self.items():
            if got >= max_njobs: break

            for rec in records:
                for accuracy in rec.ACCURACIES:
                    data = rec[accuracy]
                    # TODO: Better treatment of failed!
                    if data in  ("scheduled", "failed"): continue
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

    @add_fig_kwargs
    def plot_errors(self, reference="ae", accuracy="normal", ax=None, **kwargs):
        """
        Plot the error wrt reference

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

        #ax.scatter([high_hint], [1.0], s=20) #, c='b', marker='o', cmap=None, norm=None)

        return fig

    #def get_frame(self)


class RocksaltOutdb(GbrvOutdb):
    """Results for the rocksalt structures."""
    struct_type = "rocksalt"


class PeroviskiteOutdb(GbrvOutdb):
    """Results for the ABO3 structures."""
    struct_type = "ABO3"


class HalfHeuslersOutdb(GbrvOutdb):
    """Results for the half-Heuslers structures."""
    struct_type = "hH"


def check_consistency(json_path):
    """
    Check whether the set of results stored in the json file `json_path`
    is not consistent with the set of pseudopotentials.
    This usually happens when new pseudopotentials have been added to the dojo directory.

    Returns the number of records that have been added.
    """
    # TODO: To be tested.
    # Build the interfaces with the GBRV results and the set of dojo pseudos.
    odata = GbrvOutdb.from_file(json_path)
    dojo_pptable = DojoTable.from_dojodir(os.path.dirname(json_path))

    # This is gonna be slow!
    missing = defaultdict(list)

    for formula, species in odata.gbrv_formula_and_species:
        comb_list = dojo_pptable.all_combinations_for_elements(set(species))
        records = odata[formula]

        for pseudos in comb_list:
            for rec in records:
                if rec.matches_pseudos(pseudos): break
            else:
                missing[formula].append(pseudos)

    count = 0 
    if missing:
        for formula, pplist in missing.items():
            for pseudos in pplist:
                count += 1
                odata[formula].append(GbrvRecord(formula, pseudos, dojo_pptable))

        odata.json_write()

    return count

