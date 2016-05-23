from __future__ import unicode_literals, division, print_function

import os.path
import collections
import numpy as np
import pseudo_dojo.data as pdj_data

from pseudo_dojo.core.testing import PseudoDojoTest


class PseudoTestCase(PseudoDojoTest):

    def test_oncvpsp_dojo_report(self):
        """Testing pseudopotentials with dojo report"""
        plot = True
        try:
            from matplotlib.figure import Figure as Fig
        except ImportError:
            Fig = None
            plot = False

        h_wdr = pdj_data.pseudo("H-wdr.oncvpsp")

        # Test DOJO REPORT and md5
        assert h_wdr.symbol == "H"

        #h_wdr.check_and_fix_dojo_md5()
        ref_md5 = "0911255f47943a292c3905909f499a84"
        assert h_wdr.compute_md5() == ref_md5
        #assert "md5" in h_wdr.dojo_report and h_wdr.md5 == ref_md5

        print(repr(h_wdr))
        print(h_wdr.as_dict())

        # Test DojoReport
        report = h_wdr.read_dojo_report()
        #print(report)
        assert report.symbol == "H" and report.element.symbol == "H"
        assert not report.has_hints
        assert report["pseudo_type"] == "norm-conserving" and report["version"] == "1.0"
        assert not report.has_hints

        # Basic consistency tests.
        missings = report.find_missing_entries()
        #assert not missings
        assert "ebands" in missings
        assert "phwoa" in missings
        with self.assertRaises(report.Error): report.has_trial("foo")

        for trial in report.trials:
            assert report.has_trial(trial)
        assert report.has_trial("deltafactor", ecut=32)

        # Test deltafactor entry.
        self.assert_almost_equal(report["deltafactor"][32]["etotals"][1], -63.503524424394556)
        self.assert_almost_equal(report["deltafactor"][32]["volumes"][1],  66.80439150995784)

        #assert report.has_trial("deltafactor", ecut="32.0")
        #with self.assertRaises(report.Error): report.has_trial("deltafactor", ecut=-1)
        #with self.assertRaises(report.Error): report.has_trial("deltafactor", ecut="32.00")

        # Test GBRV entries
        self.assert_almost_equal(report["gbrv_bcc"][32]["a0"], 1.8069170394120007)
        self.assert_almost_equal(report["gbrv_fcc"][34]["a0_rel_err"], 0.044806085362549146)

        # Test Phonon entry
        self.assert_almost_equal(report["phonon"][36][-1], 528.9531110978663)

        # Test API to add ecuts and find missing entries.
        assert np.all(report.ecuts == [32.0,  34.0,  36.0, 38.0, 40.0, 42.0, 52.0])

        # TODO: reactivate these tests.
        #report.add_ecuts([30])
        #assert np.all(report.ecuts == [30.0, 32.0, 34.0, 36.0, 38.0, 40.0, 42.0, 52.0])
        #missing = report.find_missing_entries()
        #assert missing and all(v == [30] for v in missing.values())

        #report.add_ecuts([33, 53])
        #assert np.all(report.ecuts == [30.0, 32.0, 33.0, 34.0,  36.0, 38.0, 40.0, 42.0, 52.0, 53.0])
        #missing = report.find_missing_entries()
        #assert missing and all(v == [30, 33, 53] for v in missing.values())

        # Test plotting methods.
        xc = h_wdr.xc
        assert xc == "PBE"
        if plot:
            self.assertIsInstance(report.plot_deltafactor_convergence(xc=xc, show=False), Fig)
            self.assertIsInstance(report.plot_deltafactor_eos(show=False), Fig)
            self.assertIsInstance(report.plot_etotal_vs_ecut(show=False), Fig)
            self.assertIsInstance(report.plot_gbrv_convergence(show=False), Fig)
            self.assertIsInstance(report.plot_gbrv_eos('bcc', show=False), Fig)
            self.assertIsInstance(report.plot_gbrv_eos('fcc', show=False), Fig)
            self.assertIsInstance(report.plot_phonon_convergence(show=False), Fig)


#class PseudoTableTest(PymatgenTest):
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
