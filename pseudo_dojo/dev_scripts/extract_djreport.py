#!/usr/bin/env python
"""Script to extract the DOJO_REPORT from psp8 files and create external file."""
import sys
import os
import argparse
import json

from monty.os.path import find_exts
from abipy.flowtk.pseudos import Pseudo
from pseudo_dojo.core.pseudos import dojopseudo_from_file


def remove_dojo_report(path):
    """
    Remove the `DOJO_REPORT` section from the pseudopotential file.
    Write new file. Return dict with old_report, None if error.
    """
    # Read lines from file and insert jstring between the tags.
    with open(path, "rt") as fh:
        lines = fh.readlines()
        try:
            start = lines.index("<DOJO_REPORT>\n")
        except ValueError:
            start = -1

        if start == -1: return None
        stop = lines.index("</DOJO_REPORT>\n")

        report = json.loads("\n".join(lines[start+1:stop]))
        del lines[start:stop+1]

    # Write new file.
    with open(path, "w") as fh:
        fh.writelines(lines)

    return report


def main():
    try:
        top = sys.argv[1]
    except IndexError:
        print("Usage: extract_djreport TOPDIR")

    # Find all .psp8 files starting from top.
    paths = find_exts(top, ["psp8"], exclude_dirs="_*")
    #print(paths)

    for path in paths:

        try:
            pseudo = Pseudo.from_file(path)
        except Exception as exc:
            print(path, exc)
            raise

        if pseudo is None:
            print("Parser error in %s" % path)
            continue

        report_file = path.replace(".psp8", ".djrepo")
        if os.path.exists(report_file):
            #print("New DOJO file already exists. Ignoring", pseudo.filepath)
            continue

        print("Moving DOJOREPORT to", report_file)

        report = remove_dojo_report(pseudo.filepath)

        # Change md5 and pseudo_type
        report["md5"] = pseudo.compute_md5()
        if report["pseudo_type"] == "norm-conserving":
            report["pseudo_type"] = "NC"

        with open(report_file, "wt") as fh:
            json.dump(report, fh, indent=-1, sort_keys=True)
            #json.dump(report, fh, indent=4, sort_keys=True)


if __name__ == "__main__":
    sys.exit(main())
