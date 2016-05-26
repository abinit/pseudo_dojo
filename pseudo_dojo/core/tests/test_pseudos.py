from __future__ import division, print_function, unicode_literals

from pprint import pprint
from pseudo_dojo.core.testing import PseudoDojoTest
from pseudo_dojo.core.pseudos import *
from pseudo_dojo.pseudos import dojotable_absdir


class DojoTableTest(PseudoDojoTest):

    def test_from_dojodir(self):
        """Initializing DojoTable from directory."""
        table = DojoTable.from_dojodir(dojotable_absdir("ONCVPSP-PBE-PDv0.3"))

        # Produce template file for djson file.
        table.to_djson()

        # This table contains multiple pseudos for element!
        # dojo_check_errors should detect it.
        md5dict = {p.basename: p.md5 for p in table}
        errors = table.dojo_find_errors(md5dict=md5dict, require_hints=False)
        if errors: pprint(errors)
        assert errors

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
