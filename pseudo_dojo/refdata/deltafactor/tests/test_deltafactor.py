from __future__ import print_function, division

from pseudo_dojo.core import PseudoDojoTest
from pseudo_dojo.refdata.deltafactor import *


class DeltaFactorDatabaseTest(PseudoDojoTest):

    def setUp(self):
        self.database = df_database()
        new = df_database()
        self.assertTrue(new is self.database)

    def test_get_entry(self):
        """Get deltafactor entry."""
        database = self.database
        e = database.get_entry("Si")
        self.assert_almost_equal([e.v0, e.b0_GPa, e.b1], [20.4530, 88.545, 4.31])
        e = database.get_entry("Si", code="VASP")
        self.assert_almost_equal([e.v0, e.b0_GPa, e.b1], [20.446, 88.790, 4.296])

    #def test_plot(self):
    #    """Test plot_error_of_code."""
    #    self.database.plot_error_of_code("VASP", show=False)

    def test_compute_deltaf(self):
        """Computation of the delta factor."""
        database = self.database

        references = {
            "Si": 0.13551230521883104,
            "Ag": 4.1123659800975743,
        }

        for symbol, ref_df in references.items():
            print("symbol: %s", symbol)
            wien2k = database.get_entry(symbol)
            vasp = database.get_entry(symbol, code="VASP")

            df = df_compute(wien2k.v0, wien2k.b0_GPa, wien2k.b1, vasp.v0, vasp.b0_GPa, vasp.b1, b0_GPa=True)
            self.assert_almost_equal(df, ref_df, decimal=2)

            same_df = df_compute(wien2k.v0, wien2k.b0, wien2k.b1, vasp.v0, vasp.b0, vasp.b1, b0_GPa=False)
            self.assert_almost_equal(df, same_df, decimal=2)
