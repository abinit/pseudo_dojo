"""Tests for deltafactor module."""
from __future__ import print_function, division, unicode_literals

import unittest

from pymatgen.core.xcfunc import XcFunc
from pseudo_dojo.core import PseudoDojoTest
from pseudo_dojo.refdata.deltafactor.database import df_database, df_compute, read_tables_from_file


class DeltaFactorDatabaseTest(PseudoDojoTest):

    def setUp(self):
        self.pbe_db = df_database(XcFunc.from_name("PBE"))
        assert self.pbe_db.xc == "PBE"
        # Cached?
        assert df_database(XcFunc.from_abinit_ixc(11)) is self.pbe_db
        assert "WIEN2k" in self.pbe_db.codes
        assert "VASP" in self.pbe_db.codes

        spinat, spin_mode = self.pbe_db.spinat_spinmode_for_symbol("Si")
        assert spinat is None and spin_mode == "unpolarized"
        spinat, spin_mode = self.pbe_db.spinat_spinmode_for_symbol("Fe")
        assert spinat == 2 * [(0, 0, 2.3)] and spin_mode == "polarized"

        self.pw_db = df_database("PW")
        assert self.pw_db.xc == "PW"
        assert self.pw_db is not self.pbe_db
        assert "WIEN2k" in self.pw_db.codes

        # Accept PW_MOD as well.
        self.pwmod_db = df_database("PW_MOD")
        assert self.pwmod_db.xc == "PW_MOD"
        assert self.pwmod_db is not self.pw_db

    def test_get_entry(self):
        """Get deltafactor entry."""
        # PBE
        assert "Si" in self.pbe_db.symbols
        assert self.pbe_db.has_symbol("Si")
        assert not self.pbe_db.has_symbol("Foo")

        e = self.pbe_db.get_entry("Si")
        assert e.xc == self.pbe_db.xc == "PBE"
        self.assert_almost_equal([e.v0, e.b0_GPa, e.b1], [20.4530, 88.545, 4.31])
        e = self.pbe_db.get_entry("Si", code="VASP")
        self.assert_almost_equal([e.v0, e.b0_GPa, e.b1], [20.446, 88.790, 4.296])

        with self.assertRaises(self.pbe_db.Error):
            self.pbe_db.get_entry("Dy")

        # LDA-PW
        e = self.pw_db.get_entry("Si")
        assert e.xc == self.pw_db.xc == "PW"
        self.assert_almost_equal([e.v0, e.b0_GPa, e.b1], [19.694165, 96.149, 4.295])
        assert self.pw_db.get_cif_path("Si") != self.pbe_db.get_cif_path("Si")

        # PW_MOD has the same results as PW
        e_pwmod = self.pwmod_db.get_entry("Si")
        self.assert_equal([e.v0, e.b0_GPa, e.b1],
                          [e_pwmod.v0, e_pwmod.b0_GPa, e_pwmod.b1])
        assert e_pwmod.xc == "PW_MOD"
        assert e_pwmod.xc != e.xc

    def test_compute_deltaf(self):
        """Computation of the delta factor."""
        # PBE results.
        references = {
            "Si": 0.13551230521883104,
            "Ag": 4.1123659800975743,
        }

        for symbol, ref_df in references.items():
            wien2k = self.pbe_db.get_entry(symbol)
            vasp = self.pbe_db.get_entry(symbol, code="VASP")
            assert wien2k.xc == vasp.xc == "PBE"== self.pbe_db.xc

            df = df_compute(wien2k.v0, wien2k.b0_GPa, wien2k.b1, vasp.v0, vasp.b0_GPa, vasp.b1, b0_GPa=True)
            self.assert_almost_equal(df, ref_df, decimal=2)

            same_df = df_compute(wien2k.v0, wien2k.b0, wien2k.b1, vasp.v0, vasp.b0, vasp.b1, b0_GPa=False)
            self.assert_almost_equal(df, same_df, decimal=2)

    def test_plot(self):
        """Test plot_error_of_code."""
        if self.has_matplotlib():
            assert self.pbe_db.plot_error_of_code("VASP", show=False)

    #def test_read_tables_from_file(filepath):
    #    for xc, files in DeltaFactorDatabase._FILES4XC.items():
    #        for f in files:
    #            frame = read_tables_from_
