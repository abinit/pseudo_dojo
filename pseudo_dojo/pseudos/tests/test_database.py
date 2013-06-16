#!/usr/bin/env python
from __future__ import division, print_function

import os.path
import unittest

from pseudo_dojo.core import *
from pseudo_dojo.pseudos.database import pseudodojo_database, XC_TYPES, compare_checksums

##########################################################################################


class PseudoDojoDatabaseTest(PseudoDojoTest):

    def setUp(self):
        """Initialize the database."""
        self.pp_db = pseudodojo_database()

    def test_base(self):
        """Basic tests for PseudoDojoDatabase."""
        pp_db = self.pp_db

        nc_pseudos = pp_db.nc_findall_pseudos()
        self.assertTrue(nc_pseudos.allnc)

        for xc_type in XC_TYPES:
            nc_tables = pp_db.nc_tables(xc_type)
            for table in nc_tables:
                print(table)
                self.assertTrue(table.xc_type == xc_type)

        HGHK = pp_db.PBE_HGHK_TABLE
        self.assertTrue(HGHK.name == "HGHK")
        self.assertTrue(len(HGHK) == 86)
        self.assertTrue(len(HGHK.Si) == 1)

##########################################################################################

if __name__ == '__main__':
    unittest.main()
