#!/usr/bin/env python
"""Script used to rename djrepo files and associated files."""
from __future__ import division, print_function, unicode_literals

import sys
import os
import argparse


from monty.functools import prof_main
from monty.termcolor import cprint
from pseudo_dojo.core.dojoreport import DojoReport


@prof_main
def main():
    def str_examples():
        return """\
Example usage:
    djrepo_move.py src.djrepo dest.djrepo  => Move src to dst
"""

    def show_examples_and_exit(err_msg=None, error_code=0):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg:
            sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=str_examples(), formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
                        help='verbose, can be supplied multiple times to increase verbosity')
    parser.add_argument('src', nargs=1, help="Source")
    parser.add_argument('dst', nargs=1, help="Destination")

    # Parse command line.
    try:
        options = parser.parse_args()
    except Exception as exc:
        show_examples_and_exit(error_code=1)

    #if options.verbose:
    #    print("Got %s djrepo files" % len(options.paths))
    #    if options.verbose > 1:
    #        for i, p in enumerate(options.paths): print("[%d] %s" % (i, p))

    src, dst = os.path.abspath(options.src), os.path.abspath(options.dst)

    # Look before you leap.
    if not src.endswith(".djrepo") or not dst.endswith(".djrepo"):
        raise ValueError("Invalid names")
    if os.path.exists(dst):
        raise ValueError("Cannot overwrite existent djrepo file: %s" % dst)

    # Read dojoreport and change basename
    src_report = DojoReport.from_file(src)
    src_root, ext = os.path.splitext(src_report["basename"]
    src_report["basename"] = os.path.basename(dst)

    # Write new djrepo file and `git mv`
    src_report.json_write()

    #from subprocess import check_call
    def check_call(s): print("Will do ", s)
    check_call("git mv %s %s" % (src, dest))

    # `git mv` the other files.
    src_dir = os.dirname(src)
    dest_dir = os.dirname(dst)
    dest_base, _ = os.path.splitext(src)
    for f in os.listdir(src_dir):
        root, ext = os.path.splitext(f)
        if root != src_root: continue
        new = dest_base + ext
        check_call("git mv %s %s" % (os.path.join(src_dir, f), os.path.join(dest_dir, new)


if __name__ == "__main__":
    sys.exit(main())
