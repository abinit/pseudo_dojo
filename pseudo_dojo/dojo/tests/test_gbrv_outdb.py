from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os

from pseudo_dojo.core.testing import PseudoDojoTest
from pseudo_dojo.core.pseudos import DojoTable
from pseudo_dojo.dojo.gbrv_outdb import GbrvOutdb #, RocksaltOutdb, GbrvRecord,
from pseudo_dojo.pseudos import dojotable_absdir


class GbrvOutdbTest(PseudoDojoTest):

    def test_rocksalt_outdb(self):
        """Testing RocksaltOutdb database and its API."""
        return
        dirpath = dojotable_absdir("ONCVPSP-PBE")

        # Test the initialization of an empty object.
        outdb = RocksaltOutdb.new_from_dojodir(dirpath)
        #outdb.dojo_dir = "dummy_dir"
        #print(outdb)
        assert outdb.struct_type == "rocksalt"

        # Check that outdb supports pickle because the works will get a reference to it.
        self.serialize_with_pickle(outdb, protocols=None, test_eq=True)

        # Dict protocol
        assert "LiF" in outdb and "LiF" in outdb.keys()
        records = outdb["LiF"]

        # Test records (dict-like objects) supporting __eq__ and __ne__
        for rec in records:
            assert rec.formula == "LiF"
            assert rec["normal"] is None and rec["high"] is None
            assert "pseudos_metadata" in rec
            assert not rec.has_data("normal")

            d = rec.as_dict()
            assert isinstance(d, dict)
            same_rec = GbrvRecord.from_dict(d, outdb.struct_type, rec.dojo_pptable)
            #print(rec)
            assert same_rec == rec

        for formula, records in outdb.items():
            # Test find_record
            for rec in records:
                same_rec = outdb.find_record(formula, rec.pseudos)
                #assert rec.matches_pseudos(same_rec.pseudos)
                assert rec == same_rec

            # All the records for the same formula should be different!
            if len(records) > 1:
                for rec1, rec2 in zip(records[:-1], records[1:]):
                    assert rec1 != rec2

        # Here I compare all the records in the database!
        all_records = []
        for records in outdb.values():
            all_records.extend(records)

        for rec1, rec2 in zip(all_records[:-1], all_records[1:]):
            assert rec1 != rec2
            assert not rec1.matches_pseudos(rec2.pseudos)

        # Test pandas dataframe
        frame = outdb.get_dataframe()
        assert frame is not None

        # Test matplotlib tools
        if self.has_matplotlib():
            outdb.plot_errors()

        # Test API to extract jobs
        jobs = outdb.find_jobs_torun(max_njobs=3)
        assert len(jobs) == 3

        # Retrieve the record from the job params and make sure
        # the entry is set to scheduled.
        for job in jobs:
            rec = outdb.find_record(job.formula, job.pseudos)
            assert rec[job.accuracy] == "scheduled"

        # Write the object in json format
        filepath = "dummy.json"
        outdb.json_write(filepath=filepath)

        # And new we re-read it from file.
        new_outdb = GbrvOutdb.from_file(filepath)
        assert new_outdb.struct_type == outdb.struct_type
        assert len(new_outdb) == len(outdb)

        # NB: This works because all values support __eq__
        assert new_outdb == outdb

    def test_db_update(self):
        """Testing DB update"""
        return
        dirpath = dojotable_absdir("ONCVPSP-PBE")

        # Init an empty object.
        outdb = RocksaltOutdb.new_from_dojodir(dirpath)

        # No change here
        u = outdb.check_update()
        print(u)
        assert u.nrec_added == 0 and u.nrec_removed == 0

        # Now I hack a bit the object to simulate a pseudo that has been removed
        new_table = [p for p in outdb.dojo_pptable if p.basename != "Si.psp8"]
        outdb.dojo_pptable = DojoTable.as_table(new_table)

        # TODO:
        u = outdb.check_update()
        print(u)
        assert u.nrec_added == 0 and u.nrec_removed == 0
