"""Unit tests for lanthanides database."""
from __future__ import division, print_function, unicode_literals, absolute_import

from pymatgen.core.xcfunc import XcFunc

from pseudo_dojo.core.testing import PseudoDojoTest
from pseudo_dojo.refdata.lantanides.database import raren_database


class TestLaNDatabase(PseudoDojoTest):

    def test_database(self):
        """Testing lantanides database..."""
        # Init the database.
        with self.assertRaises(ValueError):
            raren_database(xc="PW")

        db = raren_database(xc="PBE")
        assert db.xc == "PBE"
        table = db.table
        assert set(table.keys()) == set("Z exp VASP VASP@600eV VASP@800eV VASP@1000eV ref Wien2k Vlab ".split())

        # Cached?
        assert db is raren_database(xc="PBE")

        # Test basic methods

        # Get FCC entry for Silicon
        #fcc_si = db.get_fcc_entry("Si")
        #assert fcc_si.ae == 3.857 and fcc_si.gbrv_uspp == 3.853
        #assert fcc_si.struct_type == "fcc" and fcc_si.species == ["Si"]
        #sfcc = fcc_si.build_structure()
        #self.assert_almost_equal(sfcc.volume, fcc_si.ae**3 /4.)

        # Get BCC entry for H
        #bcc_h = db.get_bcc_entry("H")
        #assert bcc_h.ae == 1.806 and bcc_h.gbrv_paw == 1.807
        #assert bcc_h.struct_type == "bcc" and bcc_h.species == ["H"]
        #sbcc = bcc_h.build_structure()
        #self.assert_almost_equal(sbcc.volume, bcc_h.ae**3 /2.)

        # Hg is missing.
        #missing = db.get_bcc_entry("Hg")
        #assert missing.ae is None

        # Get Rocksalt entry for LiF
        #lif = db.get_rocksalt_entry("LiF")
        #assert lif.symbol == "LiF"  and lif.struct_type == "rocksalt"
        #assert lif.species == ["Li", "F"]
        #assert lif.ae == 4.076  and lif.pslib == 4.081
        #struct = lif.build_structure()
        #assert str(struct.to_abivars())

        # Get AB03 entry for SrLiF3
        #srlif3 = db.get_abo3_entry("SrLiF3")
        #assert srlif3.symbol == "SrLiF3" and srlif3.struct_type == "ABO3"
        #assert srlif3.ae == 3.884  and srlif3.gbrv_paw == 3.883
        #struct = srlif3.build_structure()
        #assert str(struct.to_abivars())

        # Get hH entry for SrLiF3
        #agalge = db.get_hH_entry("AgAlGe")
        #assert agalge.symbol == "AgAlGe" and agalge.struct_type == "hH"
        #assert agalge.ae == 6.224  and agalge.gbrv_paw == 6.218
        #struct = agalge.build_structure()
        #assert str(struct.to_abivars())

        #assert len(db.entries_with_symbol("Si")) == 10
        #assert not db.entries_with_symbol("Lu")
        #assert db.match_symbols(["Si", "Lu"]) is None
        #assert db.match_symbols(["Si", "O"]).formula == "SiO"
        #e = db.entries_with_symbols(["Sr", "Si", "O"])
        #assert e and len(e) == 1 and e[0].formula == "SrSiO3"
