from __future__ import print_function, division

import sys
import os

from pseudo_dojo.core.testing import PseudoDojoTest
from pseudo_dojo.dojo.gbrv_outdb import GbrvRecord, GbrvOutdb, RocksaltOutdb
from pseudo_dojo.pseudos import dojo_absdir


class GbrvOutdbTest(PseudoDojoTest):

    def test_rocksalt_outdb(self):
        """Testing RocksaltOutdb database and its API."""
        dirpath = dojo_absdir("ONCVPSP-PBE")

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

            d = rec.as_dict()
            same_rec = GbrvRecord.from_dict(d, rec.dojo_pptable)
            #print(rec)
            assert same_rec == rec

        for formula, records in outdb.items():
            # Test find_record
            for rec in records:
                same_rec = outdb.find_record(formula, rec.pseudos)
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

        #records = new_outdb["ScN"]
        #pprint(records)
        #print(records)
            #print(rec)
        #assert len(records) > 1
        #rec0, rec1 = records[:2]
        #print(rec0)
        #assert rec0 == rec0
        #assert rec0.matches_pseudos(rec0.pseudos)
        #assert rec0 != rec1
        #assert not rec0.matches_pseudos(rec1.pseudos)
        #params = outdb.find_jobparams_torun()
        #print(torun)


if __name__ == "__main__":
    import unittest
    unittest.main()
