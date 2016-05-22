"""Unit tests for PseudoDojo official tables."""
from __future__ import print_function, division, unicode_literals

from pseudo_dojo.core import PseudoDojoTest
from pseudo_dojo import OfficialTables

class DojoApiTest(PseudoDojoTest):

    def test_official_tables(self):
        """Testing PseudoDojo official tables."""
        # Istantiate PseudoDojo tables.
        dojo_tables = OfficialTables()
        # Singleton?
        assert dojo_tables is OfficialTables()
        #dojo_tables.plot_numpseudos()

        # dojo_tables is a dict: table_name --> PseudoTable
        print("All dojo tables", dojo_tables)
        print("Table names", dojo_tables.keys())
        #print("available XC functionals", dojo_tables.available_xcs)

        # Can select tables by XC and or PP type.
        print("all_nctables:", dojo_tables.all_nctables())
        print("all_pawtables with xc=PBE:", dojo_tables.all_pawtables(xc="PBE"))

        # Select one particular table to be used for calculations.
        # Main entry point for client code.
        oncv_table = dojo_tables["ONCVPSP-PBE-PDv0.2-accuracy"]
        print(oncv_table)

        #assert oncv_table.xc == "PBE"
        assert all(p.isnc for p in oncv_table)
        #assert all(p.xc == oncv_table.xc for p in oncv_table)
        #frame = oncv_table.get_dfgbrv_dataframe()
        #print(frame)
        #my_table.plot_dfgbrv_dist()
        
        # Validate tables.
        retcode = 0
        #for tab in dojo_tables.values():
        #    errors = tab.dojo_find_errors(md5dict=None, require_hints=False)
        #    retcode += len(errors)
        #    if errors: print(errors)
        if retcode:
            raise RuntimeError("dojo_find_errors returned errors. Check djson files.")

        # Programmatic interface to loop over tables:
        #for xc, pp_type in product(["NC", "PAW"], ["PBE", "PW"]):
        pp_type = "NC"; xc = "PBE"
        for tab in dojo_tables.select_tables(pp_type="NC", xc="PBE"):
            for pseudo in tab:
                if pp_type == "NC": assert pseudo.isnc and not pseudo.ispaw
                if pp_type == "PAW": assert not pseudo.isnc and pseudo.ispaw
                assert pseudo.xc == xc
                #assert pseudo.xc == xc == tab.xc

        # Programmatic interface to extract pseudos.
        #assert len(dojo_tables.select_nc_pseudos("Si", xc="PBE")) == 1
        #assert len(dojo_tables.select_nc_pseudos("Si", xc="PW")) == 0
        #assert len(dojo_tables.select_paw_pseudos("Si")) == 0