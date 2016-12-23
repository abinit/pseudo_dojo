#!/usr/bin/env python
"""
Script to create a tar with the downloadable files and html for the pseudo-dojo.org website.

create the directory structure

element
    EX_version_type
    ..
    ..
..
    ..
    ..


each of these contains
    the static html version of the notebook.
    the psp8 and upf version of the pseudo
    index.html

index.html contains a list of all EX_version_type


"""
from __future__ import unicode_literals, division, print_function, absolute_import

import sys
import os
import json

def main():
    usage = "Usage: "
    if "--help" in sys.argv or "-h" in sys.argv:
        print(usage)
        return 1
    try:
        path = sys.argv[1]
    except:
        print(usage)
        return 1

    #  walk the current tree, create the directory structrure an copy the .in, .psp8, and .djrepo files

    #  walk the new tree and create the .upf and notebook html files

    #  walk the new tree again and create the index.html files.

    return

if __name__ == "__main__":
    sys.exit(main())
