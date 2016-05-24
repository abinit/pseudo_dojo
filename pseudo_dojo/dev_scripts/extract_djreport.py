#!/usr/bin/env python
"""Script to extract the DOJO_REPORT from psp8 files and create external file."""
from __future__ import division, print_function, unicode_literals

import sys
import os
import argparse
import json

#from pymatgen.io.abinit.pseudos import Pseudo
from pseudo_dojo.core.pseudos import dojopseudo_from_file
from monty.os.path import find_exts

"""
def remove_dojo_report(self):
    Remove the `DOJO_REPORT` section from the pseudopotential file.
    # Read lines from file and insert jstring between the tags.
    with open(self.path, "r") as fh:
        lines = fh.readlines()
        try:
            start = lines.index("<DOJO_REPORT>\n")
        except ValueError:
            start = -1

        if start == -1: return

        stop = lines.index("</DOJO_REPORT>\n")
        if stop == -1: return

        #del lines[start:stop]
        del lines[start:]

    # Write new file.
    with FileLock(self.path):
        with open(self.path, "w") as fh:
            fh.writelines(lines)
"""

def main():
    try:
        top = sys.argv[1]
    except IndexError:
        print("Usage: extract_djreport TOPDIR")

    # Find all .psp8 files starting from top.
    paths = find_exts(top, ["psp8"], exclude_dirs="_*")
    print(paths)

    for path in paths:
        #pseudo = Pseudo.from_file(path)
        pseudo = dojopseudo_from_file(path)
	if pseudo is None:
            print("Parser error in %s", path)
	    continue
        if not pseudo.has_dojo_report:
            print("No DOJOREPORT in %s. Ignoring file" % pseudo.filepath)
	    continue
        report_file = path.replace(".psp8", ".djrepo")
        if os.path.exists(report_file):
            print("New DOJO file already exists. Ignoring", pseudo.filepath)
	    continue
        print("Moving DOJOREPORT to %s", report_file)
        report = pseudo.read_dojo_report()

        with open(report_file, "wt") as fh:
            json.dump(report, fh, indent=-1, sort_keys=True)
            #json.dump(report, fh, indent=4, sort_keys=True)
        pseudo.remove_dojo_report()


if __name__ == "__main__":
    sys.exit(main())
