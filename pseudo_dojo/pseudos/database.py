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

##########################################################################################

PP_TYPES = ["NC", "PAW"]
                              
XC_TYPES = ["LDA", "GGA"]
                              
XC_FLAVORS  = ["PBE", "PW91"]

class DojoTable(PseudoTable):
    """
    A DojoTable is essentially a list of pseudopotential objects. 

    attributes::
        name:
            Name of the table.
        pp_type:
            Type of pseudopotential containe in the table e.g NC, PAW.
        xc_type:
            Type of XC functional e.g. LDA, GGA, MGGA.
        xc_flavor:
            String specifying the particular parametrization e.g. PW91, PBE.
        description:
            String with basic information on the table. e.g list of references.
        keywords:
            Tags associated to the table.
    """
    def __init__(self, pseudos, table_name, pp_type, xc_type, xc_flavor, description, keywords):
        """Extends PseudoTable adding metatada on the generation of the table."""
        super(DojoTable, self).__init__(pseudos)

        assert (pp_type in PP_TYPES and 
                xc_type in XC_TYPES and 
                xc_flavor in XC_FLAVORS)

        self.name = table_name
        self.pp_type = pp_type
        self.xc_type = xc_type
        self.xc_flavor = xc_flavor
        self.description = description
        self.keywords = keywords

    @classmethod
    def from_directory(cls, dirpath, metadata=None):
        exclude_exts = [".py", ".pyc", ".ini", ".sh", ".gz", ".pl", ".txt", ".swp", ".data"]
        pseudos = PseudoParser().scan_directory(dirpath, exclude_exts=exclude_exts)

        m = metadata
        if m is None:
            # Read metadata from __init__.py.
            module_name = os.path.join(dirpath, "__init__.py")
            m = imp.load_source(module_name, module_name)

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

##########################################################################################


class _PseudoDojoDatabase(dict):
    #        "TM"
    #"LDA"   "HGH"
    #"GGA"   "HGK"
    #        "FHI"
    #        "USERS"

    #"LDA"   "ATOMPAW"
    #"GGA"   "USPP"
    #        "USERS"

    def __init__(self, top=None):
        super(_PseudoDojoDatabase, self).__init__()

        # Initialize the database with empty lists. 
        for key in PP_TYPES:
            self[key] = None
                                                                        
        for pp_type in self:
            self[pp_type] = {k: [] for k in XC_TYPES}

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
                table = DojoTable.from_directory(root, metadata=m)

                # Add the table to the list
                self[m.pp_type][m.xc_type].append(table) 

    @property
    def PBE_HGHK_TABLE(self):
        for table in self.get_tables("NC", "GGA"):
            if table.name == "HGHK":
                return table
        else:
            return None

    def nc_findall_pseudos(self):
        "Return a list with all the NC pseudopotentials in the database"
        pseudos = []
        for xc_type, tables in self["NC"].items():
            for table in tables:
                pseudos.extend([p for p in table])
        return PseudoTable(pseudos)

    def get_tables(self, pp_type, xc_type):
        """Returns all the tables with the given xc_type."""
        return self[pp_type][xc_type]

    def nc_tables(self, xc_type):
        "Iterate over the norm-conserving tables with XC type xc_type."
        return self.get_tables("NC", xc_type)

    def paw_tables(self, xc_type):
        "Iterate over the PAW tables with XC type xc_type."
        return self.get_tables("PAW", xc_type)

    def nc_pseudos(self, symbol, xc_type, table_name=None, **kwargs):
        "Return a list of :class:`Pseudo` instances."
        pseudos = []
        for table in self.nc_tables(xc_type):
            if table_name is not None and table_name != table.name: continue
            pseudos.extend(table.pseudos_with_symbol(symbol))
        return PseudoTable(pseudos)

    def show(self, verbose):
        """Print basic information on the database."""
        print(self)



# Global variable storing the database of pseudopotentials.
_OFFICIAL_DATABASE = _PseudoDojoDatabase()

def reload_pseudodojo_database():
    """Reload the Database at runtime."""
    _OFFICIAL_DATABASE = _PseudoDojoDatabase()

# Official API.
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

    def __init__(self, *args, **kwargs):
        super(ChecksumsDatabase, self).__init__(*args, **kwargs)

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
            skip_exts = ["py", "pyc", "sh", "txt"]
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
        d = {k: PseudoChecksums(**d[k]) for k in d}
        return cls(**d)
                              
    @property
    def to_dict(self):
        return {k:v.to_dict for k,v in self.items()}

    @classmethod
    def json_load(cls, filename):
        print(type(cls))
        with open(filename, "r") as fp:
            d = json.load(fp)
            return cls.from_dict(d)

    def json_dump(self, filename):
        with open(filename, "w") as fp:
            json.dump(self.to_dict, fp, indent=4)

#############################################
# API for manipulating the checksum database.
#############################################

def _replace_reference_checksums(new_db):
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

def test_checksums():
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

if __name__ == "__main__":
    assert test_checksums() == 0
    import unittest
    unittest.main()
