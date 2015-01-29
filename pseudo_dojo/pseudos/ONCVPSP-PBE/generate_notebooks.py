#!/usr/bin/env python
import sys
import os

from IPython.nbformat import current as nbf

def write_notebook(pseudopath):
    """See http://nbviewer.ipython.org/gist/fperez/9716279"""
    nb = nbf.new_notebook()

    # This notebook will simply have three cells that read print 0, print 1, etc:
    cells = [
        nbf.new_text_cell('markdown', "This is an auto-generated notebook."),
        nbf.new_code_cell("""\
from __future__ import print_function
import mpld3
mpld3.enable_notebook()
import seaborn as sns
#sns.set(style="dark", palette="Set2")
sns.set(style='ticks', palette='Set2')
"""),
        nbf.new_code_cell("""\
# Construct the pseudo object and get the DojoReport
from pymatgen.io.abinitio.pseudos import Pseudo
pseudo = Pseudo.from_file('%s')
report = pseudo.dojo_report
""" % os.path.basename(pseudopath)),

        nbf.new_code_cell("""\
# Show the input file used to generate the pseudo
input_file = pseudo.filepath.replace(".psp8", ".in")
with open(input_file, "rt") as fh:
    for line in fh: print(line, end="")
"""),

        nbf.new_code_cell("""\
# Get data from the output file
from pseudo_dojo.ppcodes.oncvpsp import OncvOutputParser, PseudoGenDataPlotter
onc_parser = OncvOutputParser(pseudo.filepath.replace(".psp8", ".out"))
# Parse the file and fuild the plotter
onc_parser.scan()
plotter = onc_parser.make_plotter()
"""),

        nbf.new_code_cell("""\
# AE/PS radial wavefunctions.
plotter.plot_radial_wfs(show=False)
"""),

        nbf.new_code_cell("""\
# Arctan of the logarithmic derivatives.
plotter.plot_atan_logders(show=False)
"""),

        nbf.new_code_cell("""\
# Convergence in G-space computed by the ONCVPSP generator.
plotter.plot_ene_vs_ecut(show=False)
"""),

        nbf.new_code_cell("""\
# Projectors.
plotter.plot_projectors(show=False)
"""),

        nbf.new_code_cell("""\
# Core/Valence/Model charge densities.
plotter.plot_densities(show=False)
"""),

        nbf.new_code_cell("""\
# Local potential and l-dependent potentials.
plotter.plot_potentials(show=False)
"""),

        nbf.new_code_cell("""\
# 1-st derivatives of vl and vloc computed via finite differences.
plotter.plot_der_potentials(show=False)
"""),

        nbf.new_code_cell("""\
# Convergence of the total energy (computed from the deltafactor runs with WIEN2K equilibrium volume)
report.plot_etotal_vs_ecut(show=False)
"""),

        nbf.new_code_cell("""\
# EOS for the different cutoff energies.
report.plot_deltafactor_eos(show=False)
"""),

        nbf.new_code_cell("""\
# Convergence of the deltafactor results
report.plot_deltafactor_convergence(what=("dfact_meV", "dfactprime_meV"), show=False)
"""),

        nbf.new_code_cell("""\
# Convergence of \\Delta v0, \\Delta b0, and \\Delta b1 (deltafactor tests)
# Here we plot the difference wrt Wien2k results.
report.plot_deltafactor_convergence(what=("-dfact_meV", "-dfactprime_meV"), show=False)
"""),

        nbf.new_code_cell("""\
# Convergence of the GBRV lattice parameters.
report.plot_gbrv_convergence(show=False)
"""),
    ]

    # Now that we have the cells, we can make a worksheet with them and add it to the notebook:
    nb['worksheets'].append(nbf.new_worksheet(cells=cells))

    # Next, we write it to a file on disk that we can then open as a new notebook.
    # Note: This should be as easy as: nbf.write(nb, fname), but the current api is a little more verbose and needs a real file-like object.
    root, ext = os.path.splitext(pseudopath)
    fname = root + '.ipynb'
    with open(fname, 'w') as f:
        nbf.write(nb, f, 'ipynb')

def main():
    top = "." if len(sys.argv) > 1 else "."
    exts=("psp8",)

    def find_paths(top, exts):
        paths = []
        for dirpath, dirnames, filenames in os.walk(top):
            if os.path.basename(dirpath).startswith("_"): continue
            dirpath = os.path.abspath(dirpath)
            for filename in filenames:
                if any(filename.endswith(ext) for ext in exts):
                    paths.append(os.path.join(dirpath, filename))
        return paths

    for path in find_paths(top, exts):
        write_notebook(path)

    # Use runipy
    #from subprocess import check_output, 
    #for path in find_paths(top, exts):
    #    try:
    #        check_output(["runipy", "-o", path])
    #    except CalledProcessError as exc:
    #        print("returncode:", exc.returcode)
    #        print("output:\n%s", exc.output)

    #return 0

if __name__ == "__main__":
    sys.exit(main())
