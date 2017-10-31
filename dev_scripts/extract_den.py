#!/usr/bin/env python
from __future__ import division, print_function, unicode_literals

import sys
import os

from pseudo_dojo.ppcodes.oncvpsp import psp8_get_densities

def main():
    if "-h" in sys.argv or "--help" in sys.argv:
        print("Usage:")
        print("\textract_den.py FILE.psp8 or")
        print("\textract_den.py directory")
        return 0

    def parse_and_generate_fc(path):
        filename = os.path.basename(path)
        dirpath = os.path.dirname(path)
        if not filename.endswith(".psp8"): return 0
        path = os.path.join(dirpath, filename)
        print("About to analyse:", path)
        fc_path = os.path.join(dirpath, filename.replace(".psp8", ".fc"))
        ae_path = os.path.join(dirpath, filename.replace(".psp8", ".AE"))
        if os.path.isfile(fc_path):
            print("Warning: Cannot overwrite existent file:", fc_path)
            return 1
        if os.path.isfile(ae_path):
            print("Warning: Cannot overwrite existent file:", ae_path)
            return 1
        #fc_file, ae_file = open(fc_path, "wt"), open(ae_path, "wt")
        fc_file, ae_file = sys.stdout, sys.stdout
        psp8_get_densities(path, fc_file=fc_file, ae_file=ae_file)
        #fc_file.close()
        #ae_file.close()

        return 0

    path = os.path.abspath(sys.argv[1])
    if os.path.isfile(path):
        retcode = parse_and_generate_fc(path)
    else:
        retcode = 0
        for dirpath, dirnames, filenames in os.walk(path):
            for f in filenames:
                retcode += parse_and_generate_fc(os.path.join(dirpath, f))

    return retcode


if __name__ == "__main__":
    sys.exit(main())