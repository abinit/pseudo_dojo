#!/usr/bin/env python
from __future__ import division, print_function

import os
import unittest

from pseudo_dojo.core import *
from pseudo_dojo.pseudos.database import pseudodojo_database, XC_TYPES, ChecksumsDatabase

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



class ChecksumsDatabaseTest(PseudoDojoTest):

    def setUp(self):
        """Initialize the checksums database."""
        self.checks_db = ChecksumsDatabase().generate()

    def test_json(self):
        """Test JSON coding-encoding."""
        from tempfile import mkstemp
        tmp_fh, tmp_filename = mkstemp()

        checks_db = self.checks_db
        checks_db.json_dump(tmp_filename)

        same_db = ChecksumsDatabase.json_load(tmp_filename)
        os.remove(tmp_filename)

        self.assertTrue(type(same_db) == type(checks_db))
        self.assertTrue(set(same_db.keys()) == set(checks_db.keys()))

        err = 0
        for k,v in same_db.items():
            if same_db[k] != checks_db[k]:
                err += 1
                print(same_db[k], checks_db[k])

        self.assertTrue(err==0)


##########################################################################################

if __name__ == '__main__':
    unittest.main()
