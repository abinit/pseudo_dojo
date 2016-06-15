"""Tools to produce jupyter notebooks """
from __future__ import print_function, division, unicode_literals

import os

from monty.os.path import which


def write_notebook(pseudopath, with_validation=False, with_eos=False, tmpfile=None):
    """
    Read a pseudopotential file and write an ipython notebook.
    By default, the notebook is created in the same directory
    as pseudopath but with the extension `ipynb` unless `tmpfile` is set to True.
    In the later case, a temporay file is created.

    Args:
        pseudopath: Path to the pseudopotential file.
        with_validation: If True an ipython widget is added at the end of the notebook
          to validate the pseudopotential.
        with_eos: True if EOS plots are wanted.

    Returns:
        The path to the ipython notebook.

    See http://nbviewer.ipython.org/gist/fperez/9716279
    """
    try:
        # Here we have a deprecation warning but the API of v4 is different!
        from nbformat import current as nbf
        #import nbformat.v3 as nbf
    except ImportError:
        from IPython.nbformat import current as nbf

    nb = nbf.new_notebook()

    cells = [
        nbf.new_heading_cell("This is an auto-generated notebook for %s" % os.path.basename(pseudopath)),
        nbf.new_code_cell("""\
from __future__ import print_function, division, unicode_literals
%matplotlib notebook"""),

        nbf.new_code_cell("""\
# Construct the pseudo object and get the DojoReport
from pseudo_dojo.core.pseudos import dojopseudo_from_file
pseudo = dojopseudo_from_file('%s')
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

        nbf.new_heading_cell("Ghosts Test"),
        nbf.new_code_cell("fig = report.plot_ebands(with_soc=False, show=False)"),

        nbf.new_heading_cell("Model core charge and form factors computed by ABINIT"),
        nbf.new_code_cell("""\
with pseudo.open_pspsfile() as psps:
    psps.plot()"""),

        nbf.new_heading_cell("Convergence of the total energy:"),
        nbf.new_code_cell("""\
# Convergence of the total energy (computed from the deltafactor runs with Wien2K equilibrium volume)
fig = report.plot_etotal_vs_ecut(show=False)"""),

        nbf.new_heading_cell("Convergence of the deltafactor results:"),
        nbf.new_code_cell("""fig = report.plot_deltafactor_convergence(xc=pseudo.xc, what=("dfact_meV", "dfactprime_meV"), show=False)"""),

        nbf.new_heading_cell("Convergence of $\Delta v_0$, $\Delta b_0$, and $\Delta b_1$ (deltafactor tests)"),
        nbf.new_code_cell("""\
# Here we plot the difference wrt Wien2k results.
fig = report.plot_deltafactor_convergence(xc=pseudo.xc, what=("-dfact_meV", "-dfactprime_meV"), show=False)"""),

        nbf.new_heading_cell("deltafactor EOS for the different cutoff energies:"),
        nbf.new_code_cell("fig = report.plot_deltafactor_eos(show=False)"),

        nbf.new_heading_cell("Convergence of the GBRV lattice parameters:"),
        nbf.new_code_cell("fig = report.plot_gbrv_convergence(show=False)"),

        nbf.new_heading_cell("Convergence of phonon frequencies at $\Gamma$:"),
        nbf.new_code_cell("fig = report.plot_phonon_convergence(show=False)"),

        #nbf.new_heading_cell("Comparison with the other pseudos in this table"),
        #nbf.new_code_cell("""\
        #from pseudo_dojo import get_pseudos
        #pseudos = get_pseudos(".")
        #if len(pseudos) > 1:
        #    pseudos.dojo_compare()"""),
    ]

    if with_eos:
        # Add EOS plots
        cells.extend([
            nbf.new_heading_cell("GBRV EOS for the FCC structure:"),
            nbf.new_code_cell("""fig = report.plot_gbrv_eos(struct_type="fcc", show=False)"""),

            nbf.new_heading_cell("GBRV EOS for the BCC structure:"),
            nbf.new_code_cell("""fig = report.plot_gbrv_eos(struct_type="bcc", show=False)"""),
        ])

    if with_validation:
        # Add validation widget.
        cells.extend([
            nbf.new_heading_cell("PseudoDojo validation:"),
            nbf.new_code_cell("report.ipw_validate()"),
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


def make_open_notebook(pseudopath, with_validation=False, with_eos=True):
    """
    Generate an ipython notebook from the pseudopotential path and
    open it in the browser. Return system exit code.

    Raise:
        RuntimeError if jupyther or ipython are not in $PATH
    """
    path = write_notebook(pseudopath, with_validation=with_validation, with_eos=with_eos, tmpfile=True)

    if which("jupyter") is not None:
        return os.system("jupyter notebook %s" % path)

    if which("ipython") is not None:
        return os.system("ipython notebook %s" % path)

    raise RuntimeError("Cannot find neither jupyther nor ipython. Install them with `pip install`")
