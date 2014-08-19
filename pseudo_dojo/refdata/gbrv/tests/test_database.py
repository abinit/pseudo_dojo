"""Unit tests for gbrv database."""
from __future__ import division, print_function

from pseudo_dojo.core.testing import PseudoDojoTest
from pseudo_dojo.refdata.gbrv.database import GbrvDatabase


class TestGbrvDatabase(PseudoDojoTest):
    def test_gbrv(self):
        """Testing GBRV database..."""
        # Init the database.
        db = GbrvDatabase()

        # Test basic methods
        self.assertTrue(db.has_symbol("Si", stype="fcc"))
        self.assertFalse(db.has_symbol("Si", stype="rocksalt"))
        self.assertTrue("KMgF3" in db.all_symbols)

        # Get FCC entry for Silicon
        fcc_si = db.get_fcc_entry("Si")
        self.assertEqual(fcc_si.ae, 3.857)
        self.assertEqual(fcc_si.gbrv_uspp, 3.853)
        self.assertEqual(fcc_si.struct_type, "fcc")
        sfcc = fcc_si.build_structure()
        self.assert_almost_equal(sfcc.volume, fcc_si.ae**3 /4.)

        # Get BCC entry for H
        bcc_h = db.get_bcc_entry("H")
        self.assertEqual(bcc_h.ae, 1.806)
        self.assertEqual(bcc_h.gbrv_paw, 1.807)
        self.assertEqual(bcc_h.struct_type, "bcc")
        sbcc = bcc_h.build_structure()
        self.assert_almost_equal(sbcc.volume, bcc_h.ae**3 /2.)

        # Hg is missing.
        missing = db.get_bcc_entry("Hg")
        self.assertTrue(missing.ae is None)

        #assert 0

if __name__ == "__main__":
    import unittest
    unittest.main()
