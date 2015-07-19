# coding: utf-8
"""
Common test support for pseudo_dojo.

This single module should provide all the common functionality for tests
in a single location, so that test scripts can just import it and work right away.
"""
from __future__ import division, print_function, unicode_literals

from pymatgen.util.testing import PymatgenTest


class PseudoDojoTest(PymatgenTest):
    """Extends PymatgenTest with PseudoDojo-specific methods """

