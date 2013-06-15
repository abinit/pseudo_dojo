#!/usr/bin/env python
from __future__ import division, print_function

import os.path
import unittest

from pseudo_dojo.core import *
from pseudo_dojo.ppcodes.ape import *

##########################################################################################


class ApeInputGeneratorTest(PPDojoTest):

    def test_base(self):
        """Basic tests for ApeInputGenerator."""
        template = os.path.join(os.path.dirname(__file__), "si.apein")
        igen = ApeInputGenerator.from_template(template)

        # Test the parser.
        self.assertTrue(igen.get_scalar("Spinmode").strip() == "unpolarized")
        self.assertTrue(igen.get_block("orbitals") == 
            ['"Ne"', "3  |  0  |  2", "3  |  1  |  2", "3  |  2  |  0", "4  |  3  |  0"])

        print(igen.varnames)
        input = igen.get_strlist()

        # Test shallow copy
        new = igen.copy()
        new.set_calculationmode("pp")
        self.assertTrue(igen.get_scalar("calculationmode").strip() == "ae + pp")

##########################################################################################

if __name__ == '__main__':
    unittest.main()
