"""Unit tests for jupyter tools."""
from __future__ import unicode_literals, division, print_function

import os
import unittest
import pseudo_dojo.data as pdj_data

from pseudo_dojo.core.testing import PseudoDojoTest
from pseudo_dojo.util.notebook import write_notebook


class DojoNotebookTestCase(PseudoDojoTest):

    def test_write_notebook(self):
        """Testing write_notebook function."""
        if not self.has_nbformat():
            raise unittest.SkipTest("nbformat is not installed == skipping all tests in test_write_notebook")

        nbpath = write_notebook(pdj_data.pseudopath("Si.GGA_PBE-JTH-paw.xml"), with_eos=True, tmpfile=True)
        assert os.path.exists(nbpath)

        nbpath = write_notebook(pdj_data.pseudopath("Si.psp8"), with_eos=False, tmpfile=True)
        assert os.path.exists(nbpath)
