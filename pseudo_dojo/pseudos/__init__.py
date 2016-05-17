"""
Functions providing access to file data.
Mainly used to build the public APIs and write unit tests.
"""
from __future__ import print_function, division, unicode_literals

import os

here = os.path.dirname(__file__)

DOJOTABLE_BASEDIRS = [
    "ONCVPSP-PBE-DEV",
    "ONCVPSP-PBE-PDv0.3",
    "ONCVPSP-PW-DEV"
]


def dojotable_absdir(basedir):
    """
    Return the absolute dirpath of the table from its basename
    """
    assert basedir in DOJOTABLE_BASEDIRS
    return os.path.join(here, basedir)


#def find_ncfiles(top):
#    """
#    Find all netcdf files starting from the top-level directory top.
#    Filenames must be unique. Directories whose start with "tmp_" are
#    excluded from the search.
#
#    Returns:
#        dictionary with mapping: basename --> absolute path.
#    """
#    SILENT = 0
#    ncfiles = {}
#    for dirpath, dirnames, filenames in os.walk(top):
#
#        if "tmp_" in dirpath:
#            continue
#
#        for basename in filenames:
#            apath = os.path.join(dirpath, basename)
#            if basename.endswith(".nc"):
#
#                if basename in ncfiles:
#                    err_msg =  "Found duplicated basename %s\n" % basename
#                    err_msg += "Stored: %s, new %s\n" % (ncfiles[basename], apath)
#
#                    if not SILENT:
#                        import warnings
#                        warnings.warn(err_msg)
#                        #raise ValueError(err_msg)
#                        SILENT += 1
#
#                else:
#                    ncfiles[basename] = apath 
#
#    return ncfiles


def write_notebook(pseudopath, with_eos=False, tmpfile=None):
    """
    Write an ipython notebook to pseudopath.
    By default, the notebook is created in the same directory
    as pseudopath but with the extension `ipynb` unless `tmpfile` is set to True.
    In the later case, a temporay file is created.

    Args:
        pseudo: Path to the pseudopotential file.
        with_eos: True if EOS plots are wanted.

    Returns:
        The path to the ipython notebook.

    See http://nbviewer.ipython.org/gist/fperez/9716279
    """
    from IPython.nbformat import current as nbf
    #from IPython.nbformat import v3 as nbf
    #import IPython.nbformat as nbf

    nb = nbf.new_notebook()
    
    cells = [

        nbf.new_heading_cell("This is an auto-generated notebook for %s" % os.path.basename(pseudopath)),
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
report = pseudo.dojo_report""" % os.path.abspath(pseudopath)),

        nbf.new_heading_cell("ONCVPSP Input File:"),
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

        nbf.new_heading_cell("AE and PS radial wavefunctions $\phi(r)$:"),
        nbf.new_code_cell("fig = plotter.plot_radial_wfs(show=False)"),

        nbf.new_heading_cell("Arctan of the logarithmic derivatives:"),
        nbf.new_code_cell("fig = plotter.plot_atan_logders(show=False)"),

        nbf.new_heading_cell("Convergence in $G$-space estimated by ONCVPSP:"),
        nbf.new_code_cell("fig = plotter.plot_ene_vs_ecut(show=False)"),

        nbf.new_heading_cell("Projectors:"),
        nbf.new_code_cell("fig = plotter.plot_projectors(show=False)"),

        nbf.new_heading_cell("Core-Valence-Model charge densities:"),
        nbf.new_code_cell("fig = plotter.plot_densities(show=False)"),

        nbf.new_heading_cell("Local potential and $l$-dependent potentials:"),
        nbf.new_code_cell("fig = plotter.plot_potentials(show=False)"),

        #nbf.new_heading_cell("1-st order derivative of $v_l$ and $v_{loc}$ computed via finite differences:"),
        #nbf.new_code_cell("""fig = plotter.plot_der_potentials(order=1, show=False)"""),
        #nbf.new_heading_cell("2-nd order derivative of $v_l$ and $v_{loc}$ computed via finite differences:"),
        #nbf.new_code_cell("""fig = plotter.plot_der_potentials(order=2, show=False)"""),

        nbf.new_heading_cell("Model core charge and form factors computed by ABINIT"),
        nbf.new_code_cell("""\
with pseudo.open_pspsfile() as psps:
    psps.plot()"""),

        nbf.new_heading_cell("Convergence of the total energy:"),
        nbf.new_code_cell("""\
# Convergence of the total energy (computed from the deltafactor runs with Wien2K equilibrium volume)
fig = report.plot_etotal_vs_ecut(show=False)"""),

        nbf.new_heading_cell("Convergence of the deltafactor results:"),
        nbf.new_code_cell("""fig = report.plot_deltafactor_convergence(what=("dfact_meV", "dfactprime_meV"), show=False)"""),

        nbf.new_heading_cell("Convergence of $\Delta v_0$, $\Delta b_0$, and $\Delta b_1$ (deltafactor tests)"),
        nbf.new_code_cell("""\
# Here we plot the difference wrt Wien2k results.
fig = report.plot_deltafactor_convergence(what=("-dfact_meV", "-dfactprime_meV"), show=False)"""),

        nbf.new_heading_cell("deltafactor EOS for the different cutoff energies:"),
        nbf.new_code_cell("fig = report.plot_deltafactor_eos(show=False)"),

        nbf.new_heading_cell("Convergence of the GBRV lattice parameters:"),
        nbf.new_code_cell("fig = report.plot_gbrv_convergence(show=False)"),

        nbf.new_heading_cell("Convergence of phonon frequencies at $\Gamma$:"),
        nbf.new_code_cell("fig = report.plot_phonon_convergence(show=False)"),

#        nbf.new_heading_cell("Comparison with the other pseudos in this table"),
#        nbf.new_code_cell("""\
#from pseudo_dojo import get_pseudos
#pseudos = get_pseudos(".")
#if len(pseudos) > 1:
#    pseudos.dojo_compare()"""),

    ]

    if with_eos:
        # Add EOS plots
        cells.update([
            nbf.new_heading_cell("GBRV EOS for the FCC structure:"),
            nbf.new_code_cell("""fig = report.plot_gbrv_eos(struct_type="fcc", show=False)"""),

            nbf.new_heading_cell("GBRV EOS for the BCC structure:"),
            nbf.new_code_cell("""fig = report.plot_gbrv_eos(struct_type="bcc", show=False)"""),
        ])

    # Now that we have the cells, we can make a worksheet with them and add it to the notebook:
    nb['worksheets'].append(nbf.new_worksheet(cells=cells))

    # Next, we write it to a file on disk that we can then open as a new notebook.
    # Note: This should be as easy as: nbf.write(nb, fname), but the current api is 
    # a little more verbose and needs a real file-like object.
    if tmpfile is None:
        root, ext = os.path.splitext(pseudopath)
        nbpath = root + '.ipynb'
    else:
        import tempfile
        _, nbpath = tempfile.mkstemp(suffix='.ipynb', text=True)

    with open(nbpath, 'wt') as f:
        nbf.write(nb, f, 'ipynb')

    return nbpath


def make_open_notebook(pseudopath, with_eos=True):
    """
    Generate an ipython notebook from the pseudopotential path and  
    open it in the browser.
    """
    path = write_notebook(pseudopath, tmpfile=True)
    os.system("ipython notebook %s" % path)

