"""Unit tests for PseudoDojo official tables and public API."""
import os

from itertools import product
from collections import defaultdict
from monty.os.path import find_exts
from pymatgen.core.periodic_table import Element
from pseudo_dojo.core import PseudoDojoTest
from pseudo_dojo import OfficialTables
from pseudo_dojo.pseudos import all_dojotable_absdirs, check_pseudo_path, dojotable_absdir


# List of known extensions used for PP files.
KNOWN_PPEXTENSIONS = (".psp8", ".xml")


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
        repr(oncv_pbe_table); str(oncv_pbe_table)
        assert oncv_pbe_table
        assert oncv_pbe_table.xc == "PBE"
        assert all(p.isnc for p in oncv_pbe_table)
        assert all(p.xc == oncv_pbe_table.xc for p in oncv_pbe_table)

        # TODO: Validate tables.
        #retcode = 0
        #for tab in dojo_tables.values():
        #    errors = tab.dojo_find_errors(md5dict=None, require_hints=False)
        #    retcode += len(errors)
        #    if errors: print(errors)
        #if retcode != 0:
        #    raise RuntimeError("dojo_find_errors returned errors. Check djson files.")

        # md5 values should be unique.
        retcode = 0
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

        if retcode != 0:
            for k, v in md5_pseudo.items():
                if not isinstance(v, list): continue
                #print("md5=", v)
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
        count = 0
        #for xc, pp_type in [("PBE", "NC")]:
        for pp_type, xc in product(("NC", "PAW"), dojo_tables.available_xcs):
            for tab in dojo_tables.select_tables(pp_type=pp_type, xc=xc):
                for pseudo in tab:
                    count += 1
                    if pp_type == "NC": assert pseudo.isnc and not pseudo.ispaw
                    if pp_type == "PAW": assert not pseudo.isnc and pseudo.ispaw
                    assert pseudo.xc == xc == tab.xc
        assert count > 0

        # Programmatic interface to extract pseudos.
        assert (all(p.isnc and p.symbol == "Si" and p.xc == "PBE" for p in
                dojo_tables.select_nc_pseudos("Si", xc="PBE")))

        assert (all(p.ispaw and p.symbol == "Al" and p.xc == "PW" for p in
                dojo_tables.select_paw_pseudos("Al", xc="PW")))

        # Plots
        #if self.has_matplotlib():
            #dojo_tables.plot_numpseudos()


class DojoPseudoFilesTest(PseudoDojoTest):
    """
    Test that all the pseudopotential files and the corresponding
    djrepo files contained in the `stable` directories are consistent.
    """
    def test_all_pseudo_files(self):
        """Testing all pseudopotential files located in the stable directories."""

        def ppath_from_djpath(djpath):
            """Return the path of the pseudo from djrepo path."""
            root, ext = os.path.splitext(djpath)
            found = []
            for ext in KNOWN_PPEXTENSIONS:
                if os.path.exists(root + ext): found.append(root + ext)
            if not found:
                raise ValueError("Cannot find pseudo file associated to %s" % djpath)
            if len(found) > 1:
                raise ValueError("Found multiple files associated to %s" % djpath)
            return found[0]

        retcode, count = 0, 0
        for top in all_dojotable_absdirs():
            # Find all djrepo files
            djrepo_paths = find_exts(top, exts=("djrepo"), exclude_dirs="_*")
            #assert djrepo_paths # TODO

            for djp in djrepo_paths:
                count += 1
                try:
                    pp_path = ppath_from_djpath(djp)
                except ValueError as exc:
                    print("Exception %s\n while trying to find pp_file from djrepo: %s" % (exc, djp))
                    retcode += 1
                    continue

                try:
                    retcode += check_pseudo_path(pp_path, verbose=1)
                except Exception as exc:
                    print("Exception %s\n while testing: %s" % (exc, pp_path))
                    retcode += 1

        print("Found %s/%s (errors/totnum_files)" % (retcode, count))
        assert count > 0
        # TODO To be activated
        #assert retcode == 0

    def test_stale_files(self):
        """Seeking for stale files in the stable directories"""
        # Each item in ext_groups lists the file extensions that are supposed to be grouped together
        # For psp8 files, for example, we have the the input (.in), the output (.out)
        # the actual pseudopotential file (.psp8) and the djrepo file.
        # For PAW-XML: (.xml, .djrepo)
        # Tuples are sorted in order to speed up the search so that we don't have to test every possible permutation.
        ext_groups = set([
            tuple(sorted((".psp8", ".djrepo", ".out", ".in"))),
            tuple(sorted((".xml", ".djrepo"))),
        ])

        def find_stale_files(bnames):
            """
            This is the function that performs most of the work.
            It receives a list of basenames, constructs a one-to-many map: ext_tuple --> basenames
            removes the keywords in ext_groups, takes into account possible exceptions
            and returns the list of stale files.
            """
            # Construct mapping basename --> [ext1, ext2, ext3]
            d = defaultdict(list)
            for root, ext in (os.path.splitext(b) for b in bnames):
                d[root].append(ext)

            # Invert the mapping (one --> many)
            # (ext1, ext2, ext3)  --> list_of_basenames
            # where the the tuple on the left is sorted so that we can compare with ext_groups
            ext2basenames = defaultdict(list)
            for k, v in d.items():
                ext2basenames[tuple(sorted(v))].append(k)

            #print("before ppp ext2basenames", ext2basenames)
            out_list = []
            for k, v in ext2basenames.items():
                if k in ext_groups: continue
                out_list.extend(v)

            # Here we allow for exceptions. For example, FR pseudos might have been generated
            # with the same input file as the scalar relativistic version.
            # upf files are also ignored.
            out_list = [o for o in out_list if not (o.startswith("README") or
                       o.endswith("_r") or o.endswith("upf"))]

            return out_list

        # Now we loop over all the "stable" directories to find stale files.
        from warnings import warn
        retcode, count = 0, 0
        for top in all_dojotable_absdirs():
            if "ONCVPSP-PBE-PDv0.2" in top: continue
            for eldir in filter(os.path.exists, (os.path.join(top, e.symbol) for e in Element)):
                count += 1
                stale_files = find_stale_files(os.listdir(eldir))
                if stale_files:
                    #print("Found stale files in directory %s" % eldir)
                    warn("Found stale files in directory %s" % eldir)
                    #for i, f in enumerate(stale_files): print("\t[%d] %s" % (i, f))
                    retcode += 1

        assert count > 0
        # TODO: To be activated
        #assert retcode == 0
