#!/usr/bin/env python
"""
Script to create html version of the notebook for a give pseudo for the pseudo-dojo.org website.
"""
import sys

from pseudo_dojo.util.notebook import write_notebook_html


def main():
    pspath = sys.argv[1].replace('.in', '.psp8')
    write_notebook_html(pspath, tmpfile=False, mock=False)


if __name__ == "__main__":
    sys.exit(main())
