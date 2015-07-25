#!/usr/bin/env python
from __future__ import division, print_function, unicode_literals

import sys
import os
import argparse

from pymatgen.io.abinitio.pseudos import write_notebook

#from IPython.nbformat import current as nbf
#from IPython.nbformat import v3 as nbf

def main():
    parser = argparse.ArgumentParser(add_help=False)

    path_parser = argparse.ArgumentParser(add_help=False)
    path_parser.add_argument('top', default=".", help="Top-level directory")

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    p_generate = subparsers.add_parser('generate', parents=[path_parser], help="Generate ipython notebooks")
    p_runipy = subparsers.add_parser('runipy', parents=[path_parser], help="Execute ipython ipython with runipy")

    options = parser.parse_args()

    from monty.os.path import find_exts

    if options.command == "generate":
        # Generate ipython notebooks.
        exts=("psp8",)
        #for path in find_paths(options.top, exts):
        for path in find_exts(options.top, exts, exclude_dirs="_*|.*"):
            write_notebook(path)

    elif options.command == "runipy":
        # Use runipy to execute and update the notebook.
        # Warning: this does not work in the sense that plots are produced!
        from subprocess import check_output, CalledProcessError
        #for path in find_paths(options.top, exts="ipynb"):
        for path in find_exts(options.top, exts="ipynb", exclude_dirs="_*|.*"):
            try:
                check_output(["runipy", "-o", path])
                #check_output(["ipython", "notebook", path])
                #print("here")
            except CalledProcessError as exc:
                print("returncode:", exc.returncode)
                print("output:\n%s", exc.output)

            #from runipy.notebook_runner import NotebookRunner
            #from IPython.nbformat.current import read
            #os.chdir(os.path.dirname(path))

            #with open(path, "r") as fh:
            #    notebook = read(fh, 'json')
            #    r = NotebookRunner(notebook)
            #    r.run_notebook()
            #    from IPython.nbformat.current import write
            #    write(r.nb, open("MyOtherNotebook.ipynb", 'w'), 'json')

    else:
        raise ValueError("Don't know how to handle command %s" % options.command)

    return 0


if __name__ == "__main__":
    sys.exit(main())
