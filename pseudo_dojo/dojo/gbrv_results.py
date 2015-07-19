# coding: utf-8
"""The GBRV results for binary and ternary compunds"""
from __future__ import division, print_function, unicode_literals

#import abc
#import sys
import os
#import six
import json
#import numpy as np

from collections import OrderedDict, MutableMapping, defaultdict
from monty.io import FileLock
from monty.collections import AttrDict, dict2namedtuple
from monty.functools import lazy_property
#from monty.pprint import pprint_table
#from pymatgen.core.units import Ha_to_eV
from pymatgen.io.abinitio.pseudos import Pseudo
#from abipy.core.structure import Structure
from pseudo_dojo.core.pseudos import DojoTable
from pseudo_dojo.refdata.gbrv.database import gbrv_database, species_from_formula


import logging
logger = logging.getLogger(__name__)


class GbrvRecord(AttrDict):
    """
    Example of entry of LiCl:
        
        record = {
            formula: "LiCl",
            pseudos_meta: {
                "Li": {basename: "Li-s-high.psp8", md5: 312}, 
                "Cl": {basename: "Cl.psp8", md5: 562}],
            }
            normal: {ecut:, a:, v0: , b0:, b1:},
            high: {ecut:, a:, v0: , b0:, b1},
        }
    """
    ACCURACIES = ("normal", "high")

    @classmethod
    def from_dict(cls, d):
        d = d.copy()
        new = cls(d.pop("formula"), d.pop("pseudos_metada"))

        for acc in cls.ACCURACIES:
            try:
                new[acc] = d.pop(acc)
            except KeyError:
                pass

        assert not d
        return new

    def __init__(self, formula, pseudos_or_dict):
        """
        Initialize the record for the chemical formula and the list of 
        pseudopotentials.
        """
        keys = ("basename", "Z_val", "l_max", "md5")

        if isinstance(pseudos_or_dict, (list, tuple)):
            def get_info(p):
                """Extract the most important info from the pseudo."""
                symbol = p.symbol
                d = p.as_dict()
                return {k: d[k] for k in keys}

            meta = {p.symbol: get_info(p) for p in pseudos_or_dict}

        else:
            meta = pseudos_or_dict
            for v in meta.values():
                assert set(v.keys()) == set(keys)

        super(GbrvRecord, self).__init__(formula=formula, pseudos_metada=meta, 
                                         normal=None, high=None)

    def __eq__(self, other):
        if other is None: return False
        if self.formula != other.formula: return False
        if len(self.pseudos) != len(other.pseudos): return False

        for esym, p in self.pseudos.items():
            try:
                o = other.pseudos[esym]
            except KeyError:
                return False

            if p["md5"] != o["md5"] or p["basename"] != o["basename"]: 
                return False

        return True

    def __ne__(self, other):
        return not self == other

    def add_results(self, accuracy, results):
        # Validate input.
        assert accuracy in self.ACCURACIES
        #assert set(data.keys()) == set([ecut:, a:, v0: , b0:, b1]
        self[accuracy] = results

    def matches_pseudos(self, pseudos):
        """Return True if the list of `pseudos` matches the one in the record.""" 
        d1 = {p.symbol: p for p in self.pseudos}
        d2 = {p.symbol: p for p in pseudos}
        return d1 == d2

    #@lazy_property
    #def pseudos(self):
    #    pseudos = []
    #    for esymb, d in self.pseudos_metadata.items():
    #        path = os.path.join(top, esymb): 
    #        pseudos.append(Pseudo.from_file(path))

    #    return pseudos

    #def get_jobparams(self, accuracy):
    #    """
    #    Return a namedtuple with the paramenters to be used to run
    #    the job for the specified accuracy.
    #    """
    #    assert accuracy in self.ACCURACIES
    #    #ecut = max(p.hint_for_accuracy(accuracy) for p in self.pseudos)        
    #    #return dict2namedtuple(formula=self.formula, struct_type=?, pseudos=self.pseudos, ecut=ecut)


class GbrvResults(MutableMapping):
    """
    This object stores the results for the GBRV tests (binary and ternary compounds.
    This is a base class.
    """
    struct_type = None

    @classmethod
    def from_dojodir(cls, dojodir):
        """
        Initialize the object from a top level directory that
        contains pseudopotentials in the PseudoDojo format.
        """
        pseudos = DojoTable.from_dojodir(dojodir)

        new = cls(dojodir)
        for formula, species in new.gbrv_formula_and_species:
            # Find all the possible combinations of pseudos compatible with this formula!
            comb_list = pseudos.all_combinations_for_elements(set(species))
            #print("Found ",len(comb_list), " combinations for %s" % formula)

            for pplist in comb_list:
                new[formula].append(GbrvRecord(formula, pplist))

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

            new = cls(dojo_dir)
            for formula, dict_list in d.items():
                new[formula] = [GbrvRecord.from_dict(d) for d in dict_list]

            return new

    def __init__(self, dojodir):
        db = gbrv_database()
        ord_keys = db.tables[self.struct_type].keys()

        self.dojodir = os.path.basename(dojodir)
        self.data = OrderedDict()
        for key in ord_keys:
            self[key] = []

    def __str__(self):
        #return str(self.data)
        return self.to_json()

    @property
    def basename(self):
        return self.struct_type + ".json"

    def to_json(self):
        d = self.data.copy()
        d["struct_type"] = self.struct_type
        return json.dumps(d, indent=4, sort_keys=False)

    def to_file(self, filepath=None):
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

        records = self[record.formula]
        try:
            i = records.index(GbrvRecord(formula, pseudos))
        except ValueError:
            return None

        return records[i]

    def has_record(self, record):
        return record in self[record.formula]

    def find_jobparams_torun(self, max_num=3):
        """
        Find entries whose results have not been yet calculated.
        """
        jobs, got = [], 0

        for formula, records in self.items():
            if got >= max_num: break

            for rec in records:
                for accuracy in rec.ACCURACIES:
                    entry = getattr(rec, accuracy)
                    if got < max_num and entry not in (None, "scheduled"):
                        got += 1
                        #jobs.append(rec.get_jobparams(accuracy)
                        #settattr(rec, accuracy, "scheduled")

        # Update the database.
        if jobs:
            self.to_file()

        return jobs

    @lazy_property
    def gbrv_formula_and_species(self):
        """
        Return a list of tuples. Each tuple contains the formula and the list 
        of chemical species in the same order as in formula.
        Example: ("LiF2", ["Li", "F", "F"])
        """
        db = gbrv_database()
        formulas = db.tables[self.struct_type].keys()

        items = []
        for formula in formulas:
            items.append((formula, species_from_formula(formula)))

        return items


class RocksaltResults(GbrvResults):
    """Results for the rocksalt structures."""
    struct_type = "rocksalt"


class PeroviskiteResults(GbrvResults):
    """Results for the ABO3 structures."""
    struct_type = "ABO3"


class HalfHeuslersResults(GbrvResults):
    """Results for the half-Heuslers structures."""
    struct_type = "hH"


def check_consistency(json_path):
    """
    Check whether the set of results stored in the json file `json_path`
    is not consistent with the set of pseudopotentials.
    This usually happens when new pseudopotentials have been added to the dojo directory.

    Returns the number of records that have been added.
    """
    # Build the interfaces with the GBRV results and the set of dojo pseudos.
    results = GbrvResults.from_file(json_path)
    pp_table = DojoTable.from_dojodir(os.path.dirname(json_path))

    # This is gonna be slow!
    missing = defaultdict(list)

    for formula, species in results.gbrv_formula_and_species:
        comb_list = pp_table.all_combinations_for_elements(set(species))
        records = results[formula]

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
                results[formula].append(GbrvRecord(formula, pseudos))

        results.to_file()

    return count

