"""
This module provides a databased for accessing the GBRV results,
Client code should use the official API gbrv_database() to access the database.

Example::

    db = gbrv_database()
    fcc_si, bcc_si = db.get_fcc_bcc("Si")
    print(fcc_si.ae, fcc_si.gbrv_uspp)
"""
from __future__ import division, print_function

import sys
import os
import json

from collections import namedtuple, OrderedDict
from pymatgen.core.units import FloatWithUnit

__all__ = [
    "GBRVDatabase",
]


def count_species(formula):
    """
    Construct a counter (OrderedDict) from a chemical formula. 
    Assume chemical elements start with a capital letter.

    >>> count_species("Sn")
    OrderedDict([('Sn', 1)])

    >>> count_species("OSn")
    OrderedDict([('O', 1), ('Sn', 1)])

    >>> count_species("SnO2")
    OrderedDict([('Sn', 1), ('O', 2)])

    >>> count_species("OSnO")
    OrderedDict([('O', 2), ('Sn', 1)])
    """
    # Find positions of chemical elements.
    count, inds = OrderedDict(), []
    for i, char in enumerate(formula):
        if char.isupper(): inds.append(i)

    if not len(inds):
        raise ValueError("Chemical elements should start with a capital letter")

    # Tokenize.
    for i, start in enumerate(inds):
        if (i+1) < len(inds): 
            stop = inds[i+1]
            symbol = formula[start:stop]
        else:
            symbol = formula[start:]

        digpos = -1
        for pos, char in enumerate(symbol):
            #print(char)
            if char.isdigit():
                digpos = pos
                break

        num = 1
        if digpos != -1:
            symbol = symbol[:digpos]
            num = int(symbol[digpos:])

        if symbol in count:
            count[symbol] += num
        else:
            count[symbol] = num

    return count

from pymatgen.core.structure import Structure

class MyStructure(Structure):

    @classmethod
    def bcc(cls, a, species, **kwargs):
        """Build a bcc crystal structure."""
        lattice = np.array([
            -1,  1,  1,
             1, -1,  1,
             1,  1, -1]) * 0.5 * a

        frac_coords = np.reshape([0, 0, 0, 0.5, 0.5, 0.5], (2,3))

        return cls(lattice, species, coords, coords_are_cartesian=False, **kwargs)
                   #validate_proximity=False, to_unit_cell=False, site_properties=None)

    @classmethod
    def fcc(cls, a, species, **kwargs):
        """Build a fcc crystal structure."""
        lattice = np.array([
             1,  1,  0,
             0,  1,  1,
             1,  0, -1]) * 0.5 * a
                                                                                        
        frac_coords = np.reshape([
           0,     0,   0, 
           0.5, 0.5, 0.5,
           0.5, 0.5, 0.5,
           0.5, 0.5, 0.5], (4,3))
                                                                                        
        return cls(lattice, species, coords, coords_are_cartesian=False, **kwargs)
                   #validate_proximity=False, to_unit_cell=False, site_properties=None)


    #@classmethod
    #def NaCl(cls, a, sites, **kwargs)
    #@classmethod
    #def ABO3(cls, a, sites, **kwargs)
    #@classmethod
    #def hH(cls, a, sites, **kwargs)


class GBRVEntry(namedtuple("GBRVEntry", "symbol ae gbrv_uspp vasp pslib gbrv_paw struct_type")):
    """
    """
    def __new__(cls, **kwargs):
        """Extends the base class adding type conversion of arguments."""
        for k, v in kwargs.items():
            if k in ["symbol", "struct_type"]: continue
            if v == "-":
                # Set missing entries to None
                v = None
            else:
                v = float(v)
                #v0 = FloatWithUnit(new_args[1], "ang^3")

            kwargs[k] = v

        return super(cls, GBRVEntry).__new__(cls, **kwargs)

    def build_structure(self, from_ref="ae")
        """
        Returns the crystalline structure associated to this entry.
        Use the lattive parameters obtained with reference from_ref.
        """
        # Get structure type and lattive parameters
        stype = self.struct_type 
        a = getattr(self, from_ref)
        sites = self.sites

        if stype ==  "fcc": 
            return Structure.fcc(a, sites)

        elif stype == "bcc":
            return Structure.bcc(a, sites)

        elif stype == "NaCl":
            raise NotImplementedError()

        elif stype == "ABO3":
            raise NotImplementedError()

        elif stype == "hH":
            raise NotImplementedError()

        raise ValueError("Don't know how to construct %s structure" % stype)

    @property
    def species_counter(self):
        """Returns a dictionary chemical_elements --> number_of_atoms."""
        return count_specie(self.symbol)

    @property
    def sites(self):
        sites = []
        for symb, mult in self.species_counter.items():
            sites += mult * [symb]

        return sites


def read_table_from_file(filename):
    """
    Reads GBRV data from file filename.
    Returns a dict of `GBRVEntry` objects indexed by element symbol or chemical formula.
    """
    # File Format:
    # 0) Dict with structure type
    # 1) Comment
    # 2) Header with column names
    # 3) rows in CSV format
    # Example:
    """
    # {"struct_type": "fcc"}                                              
    # fcc testing data,,Please see supplementary materials for details.,,,
    Symbol,AE,GBRV_USPP,VASP,PSLIB,GBRV_PAW
    H,2.283,2.284,2.283,2.284,2.284
    """

    table, count = OrderedDict(), 0

    with open(filename, "r") as fh:

        for i, line in enumerate(fh):
            line = line.strip()

            if i == 0:
                # Read first line with the dictionary.
                info = json.loads(line[1:])
            elif i == 2:
                # Handle the header and costruct the names of the columns
                header = [t.lower() for t in line.split(",")]
            elif i > 2:
                # Skip comments or empty lines.
                #if line.startswith("#") or not line: continue
                tokens = line.split(",")
                assert len(header) == len(tokens)
                symbol = tokens[0] 
                # Conversion GPa --> eV / A**3
                #tokens[2] = FloatWithUnit(tokens[2], "GPa").to("eV ang^-3") 
                d = {k:v for k,v in zip(header, tokens)}
                d.update(info)
                #print(d)
                table[symbol] = GBRVEntry(**d)


    return table


class GBRVDatabaseError(Exception):
    """Exceptions raised by the database."""


class GBRVDatabase(object):
    """
    This object stores the GBRV results and provides methods to access the data.
    """
    Error = GBRVDatabaseError

    # Tables available in the database.
    table_names = [
        "fcc",
        "bcc",
        "NaCl",
        "ABO3",
        "hH",
    ]

    def __init__(self):
        """Read data from CSV files and initialize the object."""
        # Directory containing the GBRV tables in CSV format.
        datadir = os.path.abspath(os.path.dirname(__file__))
        datadir = os.path.join(datadir, "data")

        # Build dictionary of tables.
        self.tables = tables = {}

        for name in self.table_names:
            filepath = os.path.join(datadir, name + ".csv")
            tables[name] = read_table_from_file(filepath)

    #def get_table(self, name):
    #    return self.tables[name]

    @property
    def all_symbols(self):
        """Set with all symbols present in the database."""
        try:
            return self._all_symbols

        except AttributeError:
            self._all_symbols = set()
            for d in self.tables.values():
                self._all_symbols.update(d.keys())

            return self._all_symbols

    def has_symbol(self, symbol, table=None):
        """
        True if table contains the given symbol
        If table is None, we test if symbol is in the database.
        """
        if table is None:
            return symbol in self.all_symbols
        else:
            return symbol in self.tables[table]

    def get_entry(self, symbol, table):
        """
        Return `GBRVEntry` in table with the given symbol.
        None if symbol is not present

        Args:
            symbol:
                Chemical symbol or chemical formula
            table:
                String identifying the GBRV table
        """
        return self.tables[table].get(symbol, None)

    def get_all_entries(self, symbol):
        """
        Return a list with all the entries in the database 
        with the given symbol.
        """
        entries = []
        for tab_name in self.table_names:
            e = self.get_entry(symbol, tab_name)
            if e is not None: entries.append(e)
            
        return entries

    def get_fcc_bcc(self, symbol):
        """Returns the entries in the FCC and in the BCC table."""
        fcc = self.get_entry(symbol, "fcc")
        bcc = self.get_entry(symbol, "bcc")

        return fcc, bcc


# Unit tests
import unittest

class test_gbrv(unittest.TestCase):
    def test_gbrv(self):
        """Test GBRV database."""
        db = GBRVDatabase()

        self.assertTrue(db.has_symbol("Si", table="fcc"))
        self.assertFalse(db.has_symbol("Si", table="NaCl"))
        self.assertTrue("KMgF3" in db.all_symbols)

        fcc_si, bcc_si = db.get_fcc_bcc("Si")
        self.assertEqual(fcc_si.ae, 3.857)
        self.assertEqual(fcc_si.gbrv_uspp, 3.853)
        self.assertEqual(fcc_si.struct_type, "fcc")
        self.assertEqual(bcc_si.struct_type, "bcc")


if __name__ == "__main__":
    unitest.main()

