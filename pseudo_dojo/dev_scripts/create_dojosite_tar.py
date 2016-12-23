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

index.html contains a list of all EX_version_type and for each entry a link to open the static html notebook and
links to download the psp8 and upf files.
"""

from __future__ import unicode_literals, division, print_function, absolute_import

import sys
import os
import json

PSEUDOS_TO_INCLUDE = ['ONCVPSP-PBE-PDv0.3','ONCVPSP-PBEsol-PDv0.3','ONCVPSP-PW-PDv0.3']

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

    #  create the main directory for the website file-system

    if os.path.isdir('website'):
        print('directory website already exists first move it away then rerun.')
        return
    else:
        os.makedirs('website')

    #  walk the current tree, create the directory structure and copy the .in, .psp8, and .djrepo files

    for dirName, subdirList, fileList in os.walk(path):
        if 'website' in dirName or 'archive' in dirName or 'DEV' in dirName:
            continue
        else:
            print('Found directory: %s' % dirName)
            element = os.path.split(dirName)[-1]
            if len(element) > 2:  # in this case we did not hit an element but something else maybe better is element not in elements
                continue
            if not os.path.isdir(os.path.join('website', element)):
                print('creating directory for %s' % element)
                os.makedirs(os.path.join('website', element))
            pseudo = os.path.split(os.path.split(dirName)[0])[1]
            print('pseudo:', pseudo)
            if pseudo not in PSEUDOS_TO_INCLUDE:
                continue
            if not os.path.isdir(os.path.join('website', element, pseudo)):
                os.makedirs(os.path.join('website', element, pseudo))

        for fname in fileList:
            if fname.split('.')[-1] in ['in', 'psp8', 'djrepo']:
                print('\t%s' % fname)


    #  walk the new tree and create the .upf and notebook html files

    #  walk the new tree again and create the index.html files.

    return

if __name__ == "__main__":
    sys.exit(main())
