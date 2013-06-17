from __future__ import print_function, division

from pseudo_dojo.core import PseudoDojoTest
from pseudo_dojo.refdata.deltafactor import DeltaFactorDataset, compute_deltaf


class DeltaFactorDatasetTest(PseudoDojoTest):

    def setUp(self):
        self.dataset = DeltaFactorDataset()
        new = DeltaFactorDataset()
        self.assertTrue(new is self.dataset)

    def test_get_entry(self):
        """Get deltafactor entry."""
        dataset = self.dataset
        e = dataset.get_entry("Si")
        self.assert_equal([e.v0, e.b0, e.bp], [20.55, 88.846,4.329])
        e = dataset.get_entry("Si", code="VASP")
        self.assert_equal([e.v0, e.b0, e.bp], [20.446, 88.790, 4.296])

    def test_plot(self):
        """Test plot_error_of_code."""
        self.dataset.plot_error_of_code("VASP", show=False)

    def test_compute_deltaf(self):
        dataset = self.dataset
        ref = dataset.get_entry("Si")
        vasp = dataset.get_entry("Si", code="VASP")

        df = compute_deltaf(ref.v0, ref.b0, ref.bp, vasp.v0, vasp.b0, vasp.bp)
        print(df)
        assert 0
