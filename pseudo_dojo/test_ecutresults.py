from __future__ import unicode_literals, division, print_function

#import os.path
#import collections
#import numpy as np
#import pseudo_dojo.data as pdj_data

#from copy import copy
from pseudo_dojo.core.testing import PseudoDojoTest
from pseudo_dojo.core.dojoreport import *


class DojoEcutResultsTest(PseudoDojoTest):

    def test_protocol(self):
        """Testing DojoEcutResults API."""
        df_class = DojoEcutResults.class_from_name("deltafactor")
        dfres = df_class()
        assert not dfres 
        assert not dfres.get_ecuts() and not dfres.get_pawecutdgs() and not dfres.get_values("foo")

        data_10 = {"ecut": 10, "v": 1}
        data_20 = {"ecut": 20, "v": 2}
        data_30 = {"ecut": 30, "v": 3}
        dfres.insert(data_10)
        assert dfres and len(dfres) == 1
        assert dfres.get_ecuts() == [10]
        assert dfres.get_values("v") == [1]
        assert dfres.get_data_for_ecut(10) == data_10
        with self.assertRaises(ValueError): dfres.get_data_for_ecut(20)

        dfres.insert(data_30)
        dfres.insert(data_20)
        assert dfres.get_ecuts() == [10, 20, 30]
        assert dfres.get_values("v") == [1, 2, 3]
        assert dfres.get_data_for_ecut(30) == data_30

        # If an ecut is already stored, we replace the old entry.
        dup_20 = {"ecut": 20, "v": "foo"}
        dfres.insert(dup_20)
        assert len(dfres) == 3
        assert dfres.get_data_for_ecut(20) == dup_20
        assert dfres.get_ecuts() == [10, 20, 30]
        assert dfres.get_values("v") == [1, "foo", 3]

        same_dfres = DojoEcutResults.from_dict(dfres.as_dict())
        print("same_dfres", same_dfres)
        print("dfres", dfres)
        assert same_dfres.as_dict() == dfres.as_dict()
