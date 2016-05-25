"""Unit tests for PseudoDojo official tables and public API."""
from __future__ import print_function, division, unicode_literals

from itertools import product
from pseudo_dojo.core import PseudoDojoTest
from pseudo_dojo import OfficialTables


class DojoApiTest(PseudoDojoTest):

    def test_official_tables(self):
        """Testing PseudoDojo official tables."""
        # Istantiate PseudoDojo tables.
        dojo_tables = OfficialTables()
        # Singleton?
        assert dojo_tables is OfficialTables()

        # dojo_tables is a dict: table_name --> PseudoTable
        print("All dojo tables", dojo_tables)
        print("Table names", dojo_tables.keys())
        print("available XC functionals", dojo_tables.available_xcs)
        assert "PBE" in dojo_tables.available_xcs

        # Can select tables by XC and or PP type.
        print("all_nctables:", dojo_tables.all_nctables())
        print("all_pawtables with xc=PBE:", dojo_tables.all_pawtables(xc="PBE"))

        # Select one particular table to be used for calculations.
        # Main entry point for client code.
        oncv_pbe_table = dojo_tables["ONCVPSP-PBE-PDv0.2-accuracy"]
        print(oncv_pbe_table)
        assert oncv_pbe_table
        assert oncv_pbe_table.xc == "PBE"
        assert all(p.isnc for p in oncv_pbe_table)
        assert all(p.xc == oncv_pbe_table.xc for p in oncv_pbe_table)

        frame = oncv_pbe_table.get_dojo_dataframe()
        frame = oncv_pbe_table.get_dfgbrv_dataframe()
        #print(frame)
        #my_table.plot_dfgbrv_dist()
        
        # TODO: Validate tables.
        retcode = 0
        #for tab in dojo_tables.values():
        #    errors = tab.dojo_find_errors(md5dict=None, require_hints=False)
        #    retcode += len(errors)
        #    if errors: print(errors)
        if retcode:
            raise RuntimeError("dojo_find_errors returned errors. Check djson files.")

        # md5 values should be unique.
        retcode  = 0
        md5_pseudo = {} 
        for table in dojo_tables.values():
            for p in table:
                if p.md5 not in md5_pseudo:
                    md5_pseudo[p.md5] = p 
                else:
                    # We have a problem. Save the pseudo for future reference.
                    retcode += 1
                    if not isinstance(md5_pseudo[p.md5], list):
                        md5_pseudo[p.md5] = [md5_pseudo[p.md5]]
                    md5_pseudo[p.md5].append(p)

        if retcode:
            for k, v in md5_pseudo.items():
                if not isinstance(v, list): continue
                print("md5=", v)
                for p in v: print(p) 
            raise RuntimeError("Found multiple md5 in pseudo dojo tables.")

        # Test find_pseudo from md5 (in principle one could do this for each pseudo!)
        #dojo_tables.select_nc_pseudos("Si", xc="PBE"))
        p = oncv_pbe_table.pseudo_with_symbol("H")
        same_p = dojo_tables.pseudo_from_symbol_md5(p.symbol, p.md5, "NC", p.xc)
        assert p.md5 == same_p.md5
        assert p == same_p
        assert p.filepath == same_p.filepath

        # Test whether table is complete and if all pseudos have similar SOC characteristics.
        for table in dojo_tables.values():
            #assert table.is_complete(zmax=118)
            h = table.pseudo_with_symbol("H")
            assert all(p.supports_soc == h.supports_soc for p in table)

        # Programmatic interface to loop over tables:
        retcode = 0
        #for xc, pp_type in [("PBE", "NC")]:
        for pp_type, xc in product(("NC", "PAW"), dojo_tables.available_xcs):
            for tab in dojo_tables.select_tables(pp_type=pp_type, xc=xc):
                for pseudo in tab:
                    retcode += 1
                    if pp_type == "NC": assert pseudo.isnc and not pseudo.ispaw
                    if pp_type == "PAW": assert not pseudo.isnc and pseudo.ispaw
                    assert pseudo.xc == xc
                    assert pseudo.xc == xc == tab.xc
        assert retcode

        # Programmatic interface to extract pseudos.
        assert (all(p.isnc and p.symbol == "Si" and p.xc == "PBE" for p in
                dojo_tables.select_nc_pseudos("Si", xc="PBE")))

        assert (all(p.ispaw and p.symbol == "Al" and p.xc == "PW" for p in
                dojo_tables.select_paw_pseudos("Al", xc="PW")))

        # Plots
        #dojo_tables.plot_numpseudos()
