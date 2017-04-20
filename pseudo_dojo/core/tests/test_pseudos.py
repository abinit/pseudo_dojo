# coding: utf-8
from __future__ import division, print_function, unicode_literals, absolute_import

import tempfile

from pprint import pprint
from pseudo_dojo.core.testing import PseudoDojoTest
from pseudo_dojo.core.pseudos import *
from pseudo_dojo.pseudos import dojotable_absdir


class DojoTableTest(PseudoDojoTest):

    def test_from_dojodir(self):
        """Initializing DojoTable from directory."""
        table = DojoTable.from_dojodir(dojotable_absdir("ONCVPSP-PBE-PDv0.3"))
        repr(table); str(table)

        # Produce template file for djson file.
        assert isinstance(table.to_djson(), dict)

        # This table contains multiple pseudos for element!
        # dojo_check_errors should detect it.
        md5dict = {p.basename: p.md5 for p in table}
        errors = table.dojo_find_errors(md5dict=md5dict, require_hints=False)
        assert errors

        # Test Dojo DataFrame
        dojo_frame, errors = table.get_dojo_dataframe()
        #print(dojo_frame)
        # TODO
        #if errors:
        #    print("Found errors in dojotable:")
        #    pprint(errors)
        #assert not errors

        # Write ipython notebook
        if self.has_nbformat():
            table.write_notebook()

        # Test helper functions
        #dojo_frame.tabulate()

        # Test myfamilies and select_families
        assert isinstance(dojo_frame.select_family("alkaline"), dojo_frame.__class__)
        myfamilies = dojo_frame.myfamilies()
        assert myfamilies
        for family in myfamilies:
            f = dojo_frame.select_family(family)
            assert len(f) and isinstance(f, dojo_frame.__class__)

        # Test myrows and select_rows
        myrows = dojo_frame.myrows()
        assert myrows
        for row in myrows:
            f = dojo_frame.select_rows(row)
            assert len(f) and isinstance(f, dojo_frame.__class__)

        assert isinstance(dojo_frame.select_rows([1, 3]), dojo_frame.__class__)

        # Plot tools
        if self.has_matplotlib():
            dojo_frame.plot_hist(show=False)
            dojo_frame.plot_trials(show=False)

        # Test DeltaFactor, GBRV DataFrame
        dfgbrv_frame = table.get_dfgbrv_dataframe()
        if self.has_matplotlib():
            dfgbrv_frame.plot_dfgbrv_dist(show=False)

    #def test_from_djson(self):
    #    """Initializing DojoTable from djson file."""
    #    djson_path = os.path.join(dojotable_absdir("ONCVPSP-PBE-DEV"), "test_accuracy.djson")
    #    table = DojoTable.from_djson_file(djson_path)

    #    # The table must have a dojo_info dict
    #    print(table.dojo_info)
    #    assert table.dojo_info


#class PseudoTableTest(PseudoDojoTest):
#
#    def test_methods(self):
#        """Test PseudoTable methods"""
#        table = PseudoTable(ref_files("14si.pspnc",  "14si.4.hgh", "14-Si.LDA.fhi"))
#        print(table)
#        assert len(table) == 3
#        for pseudo in table:
#            assert pseudo.isnc
#        assert table.allnc and not table.allpaw
#        assert table.zlist == [14]
#
#        # Data persistence
#        self.serialize_with_pickle(table, test_eq=False)
#
#        #d = table.as_dict()
#        #PseudoTable.from_dict(d)
#        #self.assertMSONable(table)
#
#        selected = table.select_symbols("Si")
#        assert len(selected) == len(table) and selected.__class__ is table.__class__
#
#        with self.assertRaises(ValueError):
#            table.pseudos_with_symbols("Si")
