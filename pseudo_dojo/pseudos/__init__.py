from __future__ import print_function, division

import os 
import hashlib
import json

from collections import namedtuple, OrderedDict
from pprint import pprint

_TOP = os.path.abspath(os.path.dirname(__file__))


class PseudoChecksum(namedtuple("PseudoChecksum", "relpath num_lines md5")):

    @classmethod 
    def from_file(cls, path):
        """
        Return the checksum of the pseudopotential file.

        The checksum is given by the tuple (basename, num_lines, hexmd5)
        where basename if the file name, hexmd5 is the (hex) MD5 hash,
        and num_lines is the number of lines in the file.
        """
        file = os.path.abspath(path)
        hasher = hashlib.md5()
        with open(path, "r") as fh:
            text = fh.read()
            hasher.update(text)

        return PseudoChecksum(relpath=os.path.relpath(path, start=_TOP),
                              num_lines=len(text.splitlines()), 
                              md5=hasher.hexdigest(),
                             )


class ChecksumDatabase(OrderedDict):

    def __init__(self, *args, **kwargs):
        super(ChecksumDatabase, self).__init__(*args, **kwargs)

    @classmethod
    def generate(cls):

        def skip(file):
            # Skip private files 
            if file[0] in [".", "_"]: 
                return True 
            # Skip files with extension in skip_exts.
            skip_exts = ["py", "sh", "txt"]
            if "." in file and file.split(".")[-1] in skip_exts: 
                return True
                                                       
            return False
                                                       
        new = cls()
        for root, dirs, files in os.walk(_TOP):
            if "__init__.py" not in files: 
                continue
            files = [f for f in files if not skip(f)]
            for file in files:
                check = PseudoChecksum.from_file(os.path.join(root, file))
                new[check.relpath] = check

        return new

    #@classmethod
    #def from_file(cls, file):

    #def to_dict(self):
    #    return self
    #@classmethod
    #def from_dict(cls, d):
    #def json_dump(self, path):

#def get_reference_checksums():
#    return ChecksumDatabase.generate()

def generate_new_checksums():
    return ChecksumDatabase.generate()

def validate_checksums():
    new_checks = generate_new_checksums()
    ref_checks = get_reference_checksums()

    ref_keys = ref_checks.keys()
    new_keys = new.checks.keys()



if __name__ == "__main__":
    new_checks = generate_new_checksums()
    pprint(new_checks)


