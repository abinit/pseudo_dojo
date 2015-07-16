#!/usr/bin/env python
from __future__ import division, print_function, unicode_literals

import sys
import os
import argparse

from monty.string import list_strings
from IPython.nbformat import current as nbf

def write_notebook(pseudopath):
    """See http://nbviewer.ipython.org/gist/fperez/9716279"""
    nb = nbf.new_notebook()

    cells = [
        nbf.new_text_cell('heading', "This is an auto-generated notebook for %s" % os.path.basename(pseudopath)),
        nbf.new_code_cell("""\
from __future__ import print_function
%matplotlib inline
import mpld3
from mpld3 import plugins as plugs
plugs.DEFAULT_PLUGINS = [plugs.Reset(), plugs.Zoom(), plugs.BoxZoom(), plugs.MousePosition()]
mpld3.enable_notebook()
import seaborn as sns
#sns.set(style="dark", palette="Set2")
sns.set(style='ticks', palette='Set2')"""),
        nbf.new_code_cell("""\
# Construct the pseudo object and get the DojoReport
from pymatgen.io.abinitio.pseudos import Pseudo
pseudo = Pseudo.from_file('%s')
report = pseudo.dojo_report""" % os.path.basename(pseudopath)),

        nbf.new_text_cell('heading', "ONCVPSP Input File:"),
        nbf.new_code_cell("""\
input_file = pseudo.filepath.replace(".psp8", ".in") 
%cat $input_file"""),

        nbf.new_code_cell("""\
# Get data from the output file
from pseudo_dojo.ppcodes.oncvpsp import OncvOutputParser, PseudoGenDataPlotter
onc_parser = OncvOutputParser(pseudo.filepath.replace(".psp8", ".out"))
# Parse the file and build the plotter
onc_parser.scan()
plotter = onc_parser.make_plotter()"""),

        nbf.new_text_cell('heading', "AE/PS radial wavefunctions $\phi(r)$:"),
        nbf.new_code_cell("""fig = plotter.plot_radial_wfs(show=False)"""),

        nbf.new_text_cell('heading', "Arctan of the logarithmic derivatives:"),
        nbf.new_code_cell("""fig = plotter.plot_atan_logders(show=False)"""),

        nbf.new_text_cell('heading', "Convergence in $G$-space estimated by ONCVPSP:"),
        nbf.new_code_cell("""fig = plotter.plot_ene_vs_ecut(show=False)"""),

        nbf.new_text_cell('heading', "Projectors:"),
        nbf.new_code_cell("""fig = plotter.plot_projectors(show=False)"""),

        nbf.new_text_cell('heading', "Core/Valence/Model charge densities:"),
        nbf.new_code_cell("""fig = plotter.plot_densities(show=False)"""),

        nbf.new_text_cell('heading', "Local potential and $l$-dependent potentials:"),
        nbf.new_code_cell("""fig = plotter.plot_potentials(show=False)"""),

        #nbf.new_text_cell('heading', "1-st order derivative of $v_l$ and $v_{loc}$ computed via finite differences:"),
        #nbf.new_code_cell("""fig = plotter.plot_der_potentials(order=1, show=False)"""),

        #nbf.new_text_cell('heading', "2-nd order derivative of $v_l$ and $v_{loc}$ computed via finite differences:"),
        #nbf.new_code_cell("""fig = plotter.plot_der_potentials(order=2, show=False)"""),

        nbf.new_text_cell('heading', "Model core charge and form factors computed by ABINIT"),
        nbf.new_code_cell("""\n
peudo = abilab.Pseudo.from_file(pseudos[0])
with pseudo.open_pspsfile(ecut=30) as psps:
    psps.plot(what="all")"""),

        nbf.new_text_cell('heading', "Convergence of the total energy:"),
        nbf.new_code_cell("""\
# Convergence of the total energy (computed from the deltafactor runs with Wien2K equilibrium volume)
fig = report.plot_etotal_vs_ecut(show=False)"""),

        nbf.new_text_cell('heading', "Convergence of the deltafactor results:"),
        nbf.new_code_cell("""fig = report.plot_deltafactor_convergence(what=("dfact_meV", "dfactprime_meV"), show=False)"""),

        nbf.new_text_cell('heading', "Convergence of $\Delta v_0$, $\Delta b_0$, and $\Delta b_1$ (deltafactor tests)"),
        nbf.new_code_cell("""\
# Here we plot the difference wrt Wien2k results.
fig = report.plot_deltafactor_convergence(what=("-dfact_meV", "-dfactprime_meV"), show=False)"""),

        nbf.new_text_cell('heading', "deltafactor EOS for the different cutoff energies:"),
        nbf.new_code_cell("""fig = report.plot_deltafactor_eos(show=False)"""),

        nbf.new_text_cell('heading', "Convergence of the GBRV lattice parameters:"),
        nbf.new_code_cell("""fig = report.plot_gbrv_convergence(show=False)"""),

        nbf.new_text_cell('heading', "GBRV EOS for the FCC structure:"),
        nbf.new_code_cell("""fig = report.plot_gbrv_eos(struct_type="fcc", show=False)"""),

        nbf.new_text_cell('heading', "GBRV EOS for the BCC structure:"),
        nbf.new_code_cell("""fig = report.plot_gbrv_eos(struct_type="bcc", show=False)"""),

#        nbf.new_text_cell('heading', "Comparison with the other pseudos in this table"),
#        nbf.new_code_cell("""\
#from pseudo_dojo import get_pseudos
#pseudos = get_pseudos(".")
#if len(pseudos) > 1:
#    pseudos.dojo_compare()"""),

    ]

    # Now that we have the cells, we can make a worksheet with them and add it to the notebook:
    nb['worksheets'].append(nbf.new_worksheet(cells=cells))

    # Next, we write it to a file on disk that we can then open as a new notebook.
    # Note: This should be as easy as: nbf.write(nb, fname), but the current api is a little more verbose and needs a real file-like object.
    root, ext = os.path.splitext(pseudopath)
    with open(root + '.ipynb', 'w') as f:
        nbf.write(nb, f, 'ipynb')

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
