"""
This module provides a simple database to access the pseudopotential tables.
"""
from __future__ import division, print_function

import sys
import os
import imp
import collections
import json
import hashlib
import numpy as np

from pseudo_dojo.core import PseudoParser, PseudoTable

__all__ = [
    "pseudodojo_database",
]

__author__ = "Matteo Giantomassi"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"

##########################################################################################

_TOP = os.path.abspath(os.path.dirname(__file__))

# Path to the database with hash values in JSON format
_MD5CHECKS_PATH = os.path.join(_TOP, "md5checksums.json")

# Tools and helper functions.
def nested_dict_items(nested):
    """Iterate over the items of a nested Mapping (e.g. a dictionary)."""
    for (key, value) in nested.items():
        if isinstance(value, collections.Mapping):
            for (inner_key, inner_value) in nested_dict_items(value):
                yield inner_key, inner_value
        else:
            yield key, value

def compute_stats(values, labels=None):
    min_val, max_val = 2 * [values[0]]
    min_pos, max_pos = 2 * [0]

    mean, mean2 = 0.0, 0.0
    for idx, value in enumerate(values):
        mean  += value
        mean2 += value * value

        if value > max_val:
            max_val = value
            max_pos = idx

        if value < min_val:
            min_val = value
            min_pos = idx

    mean = mean / len(values)
    stdev = np.sqrt(mean2/len(values) - mean)

    stats = collections.namedtuple("stats", "mean stdev max_val max_pos min_val min_pos")
    if labels is not None:
        max_pos = labels[max_pos]
        min_pos = labels[min_pos]

    return stats(mean, stdev, max_val, max_pos, min_val, min_pos)

##########################################################################################
#
# These are the keywords used to index the tables and to query data.
# -------------------------------------------------------------------
#
_PP_TYPES = ("NC", "PAW")
                              
_XC_TYPES = ("LDA", "GGA")
                              
_XC_FLAVORS  = ("PBE", "PW91")

class DojoTable(PseudoTable):
    """
    A DojoTable is essentially a list of pseudopotential objects with some 
    handy methods inherited from `PseudoTable`.

    attributes::
        name:
            Name of the table.
        pp_type:
            Type of pseudopotential contained in the table e.g NC, PAW.
        xc_type:
            Type of XC functional e.g. LDA, GGA, MGGA.
        xc_flavor:
            String specifying the particular parametrization of for of Exc e.g. PW91, PBE.
        description:
            String with basic information on the table. e.g list of references.
        keywords:
            Tags associated to the table.
    """
    def __init__(self, pseudos, table_name, pp_type, xc_type, xc_flavor, description, keywords):
        """Extends PseudoTable adding metatada on the generation of the table."""
        super(DojoTable, self).__init__(pseudos)

        assert (pp_type in _PP_TYPES and 
                xc_type in _XC_TYPES and 
                xc_flavor in _XC_FLAVORS)

        self.name = table_name
        self.pp_type = pp_type
        self.xc_type = xc_type
        self.xc_flavor = xc_flavor
        self.description = description
        self.keywords = keywords

    def __repr__(self):
        return "<%s at %s, name = %s>" % (self.__class__.__name__, id(self), self.name)

    @classmethod
    def from_directory(cls, dirpath, metadata=None):
        """
        Initialize the object from the pseudos in directory dirpath
        if metadata is None, metada are read from the __init__.py file in dirpath.
        """
        m = metadata
        if m is None:
            # Read metadata from __init__.py.
            module_name = os.path.join(dirpath, "__init__.py")
            m = imp.load_source(module_name, module_name)

        try:
            exclude_fnames = m.exclude_fnames
        except AttributeError:
            exclude_fnames = []

        exclude_exts = [".py", ".pyc", ".ini", ".sh", ".gz", ".pl", ".txt", ".swp", ".data"]
        pseudos = PseudoParser().scan_directory(dirpath, 
                                                exclude_exts=exclude_exts, exclude_fnames=exclude_fnames)

        return cls(pseudos, m.table_name, m.pp_type, m.xc_type, m.xc_flavor, m.description, m.keywords)

    def __str__(self):
        lines = []
        app = lines.append
        app("name: %s" % self.name)
        app("pp_type: %s" % self.pp_type)
        app("xc_type: %s" % self.xc_type)
        app("xc_flavor: %s" % self.xc_flavor)
        app("description: %s" % self.description)
        app("keywords: %s" % self.keywords)
        return "\n".join(lines)

    def compare_delta_factor(self, accuracy="normal", exclude_na=True):
        """Returns a table (list of lists) comparing the deltafactor results."""
        from pseudo_dojo.refdata.deltafactor import df_database
        dfdb = df_database()

        def fmt(float):
            return "%.3f" % float

        table = []; app = table.append
        app(["pseudo", "V0", "B0", "B1", "%V0", "%B0", "%B1", "Delta"])
        deltas, names = [], []
        for pseudo in self:
            if not pseudo.has_dojo_report or "delta_factor" not in pseudo.dojo_report:
                if not exclude_na: app([pseudo.name] + 7 * ["N/A"])
                continue

            ref = dfdb.get_entry(pseudo.symbol)
            vals = pseudo.dojo_report["delta_factor"][accuracy]
            v0, b0_GPa, b1, dfact = vals["v0"], vals["b0_GPa"], vals["b1"], vals["dfact"]
            rerr_v0  = 100 * (v0 - ref.v0) / ref.v0
            rerr_b0  = 100 * (b0_GPa - ref.b0_GPa) / ref.b0_GPa
            rerr_b1  = 100 * (b1 - ref.b1) / ref.b1
            app([pseudo.name] + map(fmt, [v0, b0_GPa, b1, rerr_v0, rerr_b0, rerr_b1, dfact]))

            deltas.append(dfact)
            names.append(pseudo.name)

        stats = compute_stats(deltas, labels=names)
        print(stats)
        return table

##########################################################################################


class _PseudoDojoDatabase(dict):
    """
    Database of the official pseudopotentials of the DOJO
    We use a nested dict [pp_type][xc_type] whose entries
    are list of DojoTable instances, e.g.

    d["NC"]["GGA"] = [PBE_HGHK, PBE_FHI ...]

    To have acces to the database use the function `pseudodojo_database`.
    """
    PP_TYPES = _PP_TYPES
    XC_TYPES = _XC_TYPES
    XC_FLAVORS = _XC_FLAVORS

    def __init__(self, top=None):
        """
        Initialize the database by traversing all the directories starting from top.
        """
        super(_PseudoDojoDatabase, self).__init__()

        # Initialize the database with empty lists. 
        for key in _PP_TYPES:
            self[key] = None
                                                                        
        for pp_type in self:
            self[pp_type] = {k: [] for k in _XC_TYPES}

        # Scan all direction starting from top, parse the pseudopotential 
        # files and create the DojoTable.
        if top is None: top = _TOP
        for root, dirs, files in os.walk(top):
            if root == top or "__init__.py" not in files: 
                continue

            # Check metadata from __init__.py.
            module_name = os.path.join(root, "__init__.py")
            m = imp.load_source(module_name, module_name)

            if hasattr(m, "table_name"):
                # Hack used to disable the table
                if getattr(m, "disabled", None) == True: continue
                table = DojoTable.from_directory(root, metadata=m)

                # Add the table to the list
                self[m.pp_type][m.xc_type].append(table) 

    def _lazy_table(self, pp_type, xc_type, table_name):
        """Lazy build of table attributes."""
        hidden = "_" + table_name
        try:
            return getattr(self, hidden)
        except AttributeError:
            setattr(self, hidden, None)
            for table in self.get_tables(pp_type, xc_type):
                if table.name == table_name:
                    setattr(self, hidden, table)
                    return table

    @property
    def GGA_PBE_HGHK(self):
        return self._lazy_table("NC", "GGA", "HGHK")

    def get_tables(self, pp_type, xc_type):
        """Returns the list of tables with the given xc_type."""
        return self[pp_type][xc_type]

    def nc_tables(self, xc_type):
        """Iterate over the norm-conserving tables with XC type xc_type."""
        return self.get_tables("NC", xc_type)

    def paw_tables(self, xc_type):
        """Iterate over the PAW tables with XC type xc_type."""
        return self.get_tables("PAW", xc_type)

    def nc_pseudos(self, symbol, xc_type, table_name=None, query=None):
        """Return a `PseudoTable of NC pseudos."""
        pseudos = []
        for table in self.nc_tables(xc_type):
            if table_name is not None and table_name != table.name: continue
            pseudos.extend(table.pseudos_with_symbol(symbol))
        return PseudoTable(pseudos)

    def show(self, pp_type, xc_type=None, verbose=0):
        """Print basic information on the database."""
        pp_type = "NC"
        for xc, table in nested_dict_items(self[pp_type]):
            if xc_type is not None and xc == xc_type: 
                print(xc, table)

    def nc_findall(self, query=None):
        """Return a `PseudoTable` with all the NC pseudopotentials in the database."""
        return self.findall("NC", query=query)

    def findall(self, pp_type, query=None):
        pseudos = []
        for xc_type, tables in self[pp_type].items():
            for table in tables:
                pseudos.extend([p for p in table])

        if query is not None:
            # Build select function from query object and apply it to pseudo.
            pseudos = [p for p in pseudos if self.select(query, p)]

        return PseudoTable(pseudos)

    @staticmethod
    def select(query, pseudo):
        """
        Example::
            {"l_max": 2}, {"zval": { ".in.": [1:3]}}
        """
        for k, v in query.items():

            if isinstance(v, dict):
                if v:
                    raise NotImplementedError("")
                else:
                    gotit = hasattr(pseudo, k)

            else:
                try:
                    gotit = (getattr(pseudo, k) == v)

                except AttributeError:
                    return False

            if not gotit:
                return False

        return gotit

########################################################
# Official API.
########################################################

# Singleton-like instance.
_OFFICIAL_DATABASE = _PseudoDojoDatabase()

def reload_pseudodojo_database():
    """Reload the Database at runtime."""
    _OFFICIAL_DATABASE = _PseudoDojoDatabase()

def pseudodojo_database():
    """Returns an instance of `PseudoDojoDatabase`."""
    return _OFFICIAL_DATABASE

##########################################################################################


#TODO this object can be moved to abinitio.
class PseudoChecksums(collections.namedtuple("PseudoChecksums", "md5 ppdata_md5 dojo_md5")):
    """
    This object stores the checksums associated to the pseudopotential file. 
    We use three different hash values:

        * md5: hash value associated to the entire file
        * ppdata_md5: hash value associated to the pseudopotential data i.e. the section
                      that precedes the <DOJO_REPORT> line
        * dojo_md5: hash value associated to the <DOJO_REPORT> section.
    """
    @classmethod 
    def from_file(cls, path):
        """
        Returns an instance of `PseudoChecksums from the path of the pseudopotential.
        """
        path = os.path.abspath(path)
        with open(path, "r") as fh:
            text = fh.read()

        hasher = hashlib.md5()
        hasher.update(text) 
        md5 = hasher.hexdigest()

        idx = text.find("<DOJO_REPORT>")
        ppdata_hasher = hashlib.md5()
        dojo_hasher = hashlib.md5()

        if idx == -1:
            # no dojo_report section => set dojo_md5 to None.
            ppdata_hasher.update(text) 
            ppdata_md5 = ppdata_hasher.hexdigest()
            dojo_md5 = None
        else:
            ppdata_hasher.update(text[:idx-1])
            ppdata_md5 = ppdata_hasher.hexdigest()
            dojo_hasher.update(text[idx:])
            dojo_md5 = dojo_hasher.hexdigest()

        return cls(md5=md5, ppdata_md5=ppdata_md5, dojo_md5=dojo_md5)

    @classmethod
    def from_dict(cls, d):
        return cls(**d)

    @property
    def to_dict(self):
        return self._asdict()



class ChecksumsDatabase(collections.OrderedDict):
    """OrderedDict with the hash values of the pseudopotentials."""

    #def __init__(self, *args, **kwargs):
    #    super(ChecksumsDatabase, self).__init__(*args, **kwargs)

    @classmethod
    def generate(cls):
        """
        Generate a new instance by computing the hash values 
        of the pseudopotential files stored in this subpackage.
        """
        def skip(file):
            # Skip private files 
            if file[0] in [".", "_"] or file == os.path.basename(_MD5CHECKS_PATH):
                return True 

            # Skip files with extension in skip_exts.
            skip_exts = ["py", "pyc", "sh", "txt", "old"]
            if "." in file and file.split(".")[-1] in skip_exts: 
                return True
                                                       
            return False
                                                       
        new = cls()
        for root, dirs, files in os.walk(_TOP):
            if "__init__.py" not in files: 
                continue
            files = [f for f in files if not skip(f)]
            for file in files:
                path = os.path.join(root, file)
                relpath=os.path.relpath(path, start=_TOP)
                checks = PseudoChecksums.from_file(path)
                new[relpath] = checks

        return new

    @classmethod
    def from_dict(cls, d):
        """Init an instance from a dictionary."""
        d = {k: PseudoChecksums(**d[k]) for k in d}
        return cls(**d)
                              
    @property
    def to_dict(self):
        """JSON reprensentation."""
        return {k:v.to_dict for k,v in self.items()}

    @classmethod
    def json_load(cls, filename):
        """Initialize an instance from the JSON data stored in filename."""
        with open(filename, "r") as fp:
            d = json.load(fp)
            return cls.from_dict(d)

    def json_dump(self, filename):
        """Save the object in JSON format."""
        with open(filename, "w") as fp:
            json.dump(self.to_dict, fp, indent=4)

#############################################
# API for manipulating the checksum database.
#############################################

def replace_reference_checksums(new_db):
    """Replace the reference hash table with new_db."""
    # Keep a backup copy.
    import shutil
    shutil.copy(_MD5CHECKS_PATH, _MD5CHECKS_PATH + ".old")
    # Dump new file.
    new_db.json_dump(_MD5CHECKS_PATH)

def get_reference_checksums():
    """Return the database of reference hash values."""
    return ChecksumsDatabase.json_load(_MD5CHECKS_PATH)

def get_new_checksums():
    """Genererate a new database of hash values from the pseudopotential files."""
    return ChecksumsDatabase.generate()

def regenerate_checksums(verbose=0):
    print("Regenerating checksums")
    new_checks = get_new_checksums()
    replace_reference_checksums(new_checks)

def compare_checksums():
    """
    Validate the reference hash table with the one generated from the pseudopotential files.

    Returns: (changed, hask_check) where
        changed is the number of files that are changed.
        hash_check is a namedtuple with the list of files that have been (removed, added, modified).
    """
    ref_checks = get_reference_checksums()
    new_checks = get_new_checksums()

    removed, added, modified = [], [], []

    for ref_rpath, ref_check in ref_checks.items():
        if ref_rpath not in new_checks:
            removed.append(ref_rpath)
        else:
            if ref_check != new_checks[ref_rpath]:
                modified.append(ref_rpath)

    for new_rpath, check in new_checks.items():
        if new_rpath not in ref_checks:
            added.append(new_rpath)

    hash_check = collections.namedtuple("HashCheck", "removed added modified")(removed, added, modified)
    changed = sum(map(len, hash_check))
    return changed, hash_check

def validate_checksums(verbose=0):
    """Validating checksum table."""
    err = 0

    changed, hash_check = compare_checksums()
    if not changed:
        return err

    if hash_check.removed:
        err = -2
        print("Warning: removed pseudos: %s" % hash_check.removed)

    if hash_check.added:
        err = -1
        print("Comment: added pseudos: %s" % hash_check.added)

    if len(hash_check.modified):
        err = 1
        print("Warning: modified pseudos %s" % hash_check.modified)

    return err

##########################################################################################
