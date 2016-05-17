#!/usr/bin/env python
from __future__ import division, print_function, unicode_literals

import collections
import os

from pseudo_dojo.core import *
from pseudo_dojo.refdata.nist import database as nist_database


class AtomicConfigurationTest(PseudoDojoTest):

    def test_neutrals(self):
        """Neutral configurations."""
        for symbol in nist_database.allsymbols:
            aconf = AtomicConfiguration.neutral_from_symbol(symbol)

            print(aconf)
            for state in aconf:
                print(state)

            # Test the shallow copy of states.
            newc = aconf.copy()
            self.assertFalse(newc.states is aconf.states) 
            self.assertTrue(newc.states == aconf.states) 
            self.assertTrue(aconf == newc)

            newc.add_state(n=15, l=0, occ=1.0)
            self.assertTrue(aconf != newc)

            # Cannot add the same state twice.
            with self.assertRaises(ValueError):
                newc.add_state(n=15, l="s", occ=1.0)

            # aconf is neutral, newc is not.
            self.assertEqual(symbol, aconf.symbol)
            self.assertTrue(aconf.isneutral)
            self.assertFalse(newc.isneutral)

            self.assertTrue(aconf.spin_mode == "unpolarized")

            # Remove the new state.
            newc.remove_state(n=15, l="s", occ=1.0)
            self.assertTrue(newc.isneutral)
            self.assertTrue(aconf == newc)

            # Cannot remove the same state twice
            with self.assertRaises(ValueError):
                newc.remove_state(n=15, l="s", occ=1.0)

    def test_initfromstring(self):
        """Initialization of atomic configurations from string"""
        for (symbol, confstr) in nist_database._neutral.items():
            print("symbol",symbol, confstr)
            Z = nist_database.Z_from_symbol(symbol)
            aconf = AtomicConfiguration.from_string(Z, confstr)

            self.assertTrue(aconf == AtomicConfiguration.neutral_from_symbol(symbol))


class RadialFunctionTest(PseudoDojoTest):

    def test_base(self):
        """Basic tests for RadialFunction."""
        filename = os.path.join(os.path.dirname(__file__), "wf-3s.ape")

        rf = RadialFunction.from_filename(filename)
        rf_der = RadialFunction.from_filename(filename, cols=(0,2)) 

        self.assertTrue(isinstance(rf, collections.Iterable))
        print(self)

        # Integral in 3D
        self.assert_almost_equal(rf.integral3d(), 1.0)

        for r, v in rf:
            rf.derivatives(r)

        rslice, vslice = rf[1:4]
        self.assert_equal(rslice, rf.rmesh[1:4])
        self.assert_equal(vslice, rf.values[1:4])


if __name__ == '__main__':
    import unittest
    unittest.main()
