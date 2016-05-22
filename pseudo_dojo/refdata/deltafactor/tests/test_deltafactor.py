from __future__ import print_function, division

from pseudo_dojo.core import PseudoDojoTest
from pseudo_dojo.refdata.deltafactor import df_database, df_compute


class DeltaFactorDatabaseTest(PseudoDojoTest):

    def setUp(self):
        self.pbe_db = df_database()
        self.pbe_db.xc == "PBE"
        # Cached?
        self.assertTrue(df_database() is self.pbe_db)

        self.pw_db = df_database("PW")
        assert self.pw_db.xc == "PW"
        assert self.pw_db is not self.pbe_db

    def test_get_entry(self):
        """Get deltafactor entry."""
        # PBE
        assert "Si" in self.pbe_db.symbols
        # FIXME
        #assert self.pbe_db.has_symbol("Si")
        #assert not self.pbe_db.has_symbol("Foo")

        e = self.pbe_db.get_entry("Si")
        assert e.xc == self.pbe_db.xc == "PBE"
        self.assert_almost_equal([e.v0, e.b0_GPa, e.b1], [20.4530, 88.545, 4.31])
        e = self.pbe_db.get_entry("Si", code="VASP")
        self.assert_almost_equal([e.v0, e.b0_GPa, e.b1], [20.446, 88.790, 4.296])

        # LDA-PW
        e = self.pw_db.get_entry("Si")
        assert e.xc == self.pw_db.xc == "PW"
        self.assert_almost_equal([e.v0, e.b0_GPa, e.b1], [19.694165, 96.149, 4.295])
        assert self.pw_db.get_cif_path("Si") != self.pbe_db.get_cif_path("Si")

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

    #def test_plot(self):
    #    """Test plot_error_of_code."""
    #    self.pbe_db.plot_error_of_code("VASP", show=False)
