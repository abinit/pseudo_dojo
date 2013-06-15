from __future__ import print_function, division

from pseudo_dojo.core import PseudoDojoTest
from pseudo_dojo.refdata.deltafactor import DeltaFactorDataset


class DeltaFactorDatasetTest(PseudoDojoTest):

    def setUp(self):
        self.dataset = DeltaFactorDataset()
        new = DeltaFactorDataset()
        self.assertTrue(new is self.dataset)

    def test_get_entry(self):
        dataset = self.dataset
        e = dataset.get_entry("Si")
        self.assert_equal([e.v0, e.b0, e.bp], [20.55, 88.846,4.329])
        e = dataset.get_entry("Si", code="VASP")
        self.assert_equal([e.v0, e.b0, e.bp], [20.446, 88.790, 4.296])

    def test_plot(self):
        self.dataset.plot_error_of_code("VASP", show=False)
