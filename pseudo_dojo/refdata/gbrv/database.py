"""
This module provides a databased for accessing the GBRV results,
Client code should use the official API gbrv_database() to access the database.

Example::

    db = gbrv_database()
    fcc_si = db.get_fcc_entry("Si")
    print(fcc_si.ae, fcc_si.gbrv_uspp)
"""
from __future__ import print_function, division, unicode_literals

import os
import json
import numpy as np

from collections import namedtuple, OrderedDict
from pymatgen.core.units import FloatWithUnit
from abipy.core.structure import Structure

__all__ = [
    "gbrv_database",
]


def count_species(formula):
    """
    Construct a counter (OrderedDict) from a chemical formula. 
    Assume chemical symbols start with a capital letter.
    The order of the symbols in formula is maintained.

    >>> count_species("Sn") == OrderedDict([('Sn', 1)])
    True
    >>> count_species("OSn") == OrderedDict([('O', 1), ('Sn', 1)])
    True
    >>> count_species("SnO2") == OrderedDict([('Sn', 1), ('O', 2)])
    True
    >>> count_species("OSnO") == OrderedDict([('O', 2), ('Sn', 1)])
    True
    """
    # Find positions of chemical elements.
    count, inds = OrderedDict(), []
    for i, char in enumerate(formula):
        if char.isupper():
            inds.append(i)

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
            if char.isdigit():
                digpos = pos
                break

        num = 1
        if digpos != -1:
            num = int(symbol[digpos:])
            symbol = symbol[:digpos]

        # Accumulate.
        if symbol in count:
            count[symbol] += num
        else:
            count[symbol] = num

    return count


class GbrvEntry(namedtuple("GbrvEntry", "symbol ae gbrv_uspp vasp pslib gbrv_paw struct_type")):
    """
    Store the GBRV lattice parameter obtained with the different codes and pseudos.
    Missing values are replaced by None

    Attributes:
        symbol: Chemical symbol of formula.
        ae:  AE results
        gbrb_uspp:  Ultra-soft PP results (Espresso code)
        vasp:  VASP PAW results
        pslib:  Espresso PAW results
        gbrv_paw:  PAW results (Abinit code)
        struct_type:  Structure type used to select the appropriate table 
            possible values are listed in GbrvDatabase.all_struct_types.

    .. note::
        Lattice parameters are in Angstrom
    """
    def __new__(cls, **kwargs):
        """Extends the base class adding type conversion of arguments."""
        for k, v in kwargs.items():
            if k in ["symbol", "struct_type"]: continue
            if v == "-":
                # Set missing entries to None
                v = None
            else:
                # Values in GBRV tables are in Angstrom.
                v = FloatWithUnit(v, "ang")

            kwargs[k] = v

        return super(cls, GbrvEntry).__new__(cls, **kwargs)

    def build_structure(self, ref="ae"):
        """
        Returns the crystalline structure associated to this entry.
        Use the optimized lattice parameters obtained from reference ref.
        Returns None if no lattice parameter is available.
        """
        # Get structure type and lattice parameters
        stype, a = self.struct_type, getattr(self, ref)

        # Handle missing value.
        if a is None:
            return None

        if stype == "bcc":
            return Structure.bcc(a, species=[self.symbol])

        elif stype == "fcc":
            return Structure.fcc(a, species=[self.symbol])

        elif stype == "rocksalt":
            raise NotImplementedError()
            #return Structure.rocksalt(a, sites)

        elif stype == "ABO3":
            raise NotImplementedError()
            #return Structure.peroviskite(a, sites)

        elif stype == "hH":
            raise NotImplementedError()
            #return Structure.hH(a, sites)

        raise ValueError("Don't know how to construct %s structure" % stype)

    @property
    def specie_counter(self):
        """Returns a dictionary chemical_elements --> number_of_atoms."""
        try:
            return self._specie_counter

        except AttributeError:
            self._specie_counter = count_species(self.symbol)
            return self._specie_counter

    @property
    def ntypat(self):
        """Number of type of atoms."""
        return len(self.specie_counter)


def read_table_from_file(filename):
    """
    Reads GBRV data from file filename.

    Returns a dict of :class:`GbrvEntry` objects indexed by element symbol or chemical formula.

    File Format:
        0) Dict with structure type
        1) Comment
        2) Header with column names
        3) rows in CSV format
    
    Example:
    
        # {"struct_type": "fcc"}                                              
        # fcc testing data,,Please see supplementary materials for details.,,,
        # Symbol,AE,GBRV_USPP,VASP,PSLIB,GBRV_PAW
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
                if line.startswith("#") or not line: continue
                tokens = line.split(",")
                assert len(header) == len(tokens)
                symbol = tokens[0] 

                d = {k: v for k, v in zip(header, tokens)}
                d.update(info)
                table[symbol] = GbrvEntry(**d)

    return table


class GbrvDatabase(object):
    """
    This object stores the GBRV results and provides methods to access the data.
    """
    # List of structure types used to index the tables available in the database.
    all_struct_types = [
        "fcc",
        "bcc",
        "rocksalt",
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

        for stype in self.all_struct_types:
            filepath = os.path.join(datadir, stype + ".csv")
            tables[stype] = read_table_from_file(filepath)

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

    def has_symbol(self, symbol, stype=None):
        """
        True if the table associated to structure type stype contains the given symbol
        If stype is None, we test if symbol is in the database.
        """
        if stype is None:
            return symbol in self.all_symbols
        else:
            return symbol in self.tables[stype]

    def get_entry(self, symbol, stype):
        """
        Return :class:`GbrvEntry` in the table associated to structure type stype 
        and with the given symbol. Return None if symbol is not present.

        Args:
            symbol: Chemical symbol or chemical formula
            stype: Structure type identifying the GBRV table
        """
        return self.tables[stype].get(symbol, None)

    def get_all_entries(self, symbol):
        """
        Return a list with all the entries in the database with the given symbol.
        """
        entries = []
        for stype in self.all_struct_types:
            e = self.get_entry(symbol, stype)
            if e is not None: entries.append(e)
            
        return entries

    def get_fcc_entry(self, symbol):
        """Returns the entry in the FCC table."""
        return self.get_entry(symbol, "fcc")

    def get_bcc_entry(self, symbol):
        """Returns the entry in the BCC table."""
        return self.get_entry(symbol, "bcc")



###################################
# Public API to access the database
###################################

__GBRV_DATABASE = None


def gbrv_database():
    """Returns the GBRV database with the reference results."""
    global __GBRV_DATABASE
    if __GBRV_DATABASE is None:
        __GBRV_DATABASE = GbrvDatabase()
    return __GBRV_DATABASE
