"""Tools to produce jupyter notebooks """
from __future__ import print_function, division, unicode_literals

import os
import io

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

    See also:
        http://nbviewer.jupyter.org/github/maxalbert/auto-exec-notebook/blob/master/how-to-programmatically-generate-and-execute-an-ipython-notebook.ipynb
    """
    import nbformat
    nbf = nbformat.v4
    nb = nbf.new_notebook()

    nb.cells.extend([
        nbf.new_markdown_cell("# PseudoDojo notebook for %s" % os.path.basename(pseudopath)),
        nbf.new_code_cell("""\
from __future__ import print_function, division, unicode_literals
%matplotlib notebook"""),

        nbf.new_markdown_cell("## Construct the pseudo object and the DojoReport"),

        nbf.new_code_cell("""\
from pseudo_dojo.core.pseudos import dojopseudo_from_file
pseudo = dojopseudo_from_file('%s')
report = pseudo.dojo_report""" % os.path.abspath(pseudopath)),

        nbf.new_markdown_cell("## ONCVPSP Input File"),
        nbf.new_code_cell("""\
input_file = pseudo.filepath.replace(".psp8", ".in")
%cat $input_file"""),

        nbf.new_code_cell("""\
# Get data from the oncvpsp output file
from pseudo_dojo.ppcodes.oncvpsp import OncvOutputParser, PseudoGenDataPlotter
onc_parser = OncvOutputParser(pseudo.filepath.replace(".psp8", ".out"))

# Parse the file and build the plotter
onc_parser.scan()
plotter = onc_parser.make_plotter()"""),

        nbf.new_markdown_cell("## AE and PS radial wavefunctions $\phi(r)$"),
        nbf.new_code_cell("fig = plotter.plot_radial_wfs(show=False)"),

        nbf.new_markdown_cell("""## Arctan of the logarithmic derivatives

From oncvpsp documentation:
The plots show $\phi(E) = \\arctan(R * d \psi(r)/dr |_R)$ for some $R$
greater than the core radius, where $\psi$ is the solution of the non-local
radial equation regular at the origin (i.e., the outward-integrated solution).
For a well-designed pseudopotential, $\phi(E)$ will closely track that of the all-electron potential
over a wide range of energies from well-below to well-above the valence semi-core states of interest.
The steps of $\pi$ indicate localized pseudo wave functions.
Spurious steps of $\pi$ indicate "ghost" states, which are localized states than on investigation
turn out to have more nodes than appropriate for their energies.

For $GW$ pseudos, no significant deviation should be present up to 8 Hartree."""),
        nbf.new_code_cell("fig = plotter.plot_atan_logders(show=False)"),

        nbf.new_markdown_cell("""## Convergence in $G$-space estimated by ONCVPSP
These results are obtained in the atomic configuration and should give a reasonable estimate
of the convergence behaviour wrt to `ecut` in crystalline systems."""),
        nbf.new_code_cell("fig = plotter.plot_ene_vs_ecut(show=False)"),

        nbf.new_markdown_cell("""## Projectors

In general the second projector in any channel should have one node more that the first one.
Pushing the energy of the second projector too high may cause an additional node.
This will most likely introduce ghosts."""),
        nbf.new_code_cell("fig = plotter.plot_projectors(show=False)"),

        nbf.new_markdown_cell("""## Core-Valence-Model charge densities

Much better convergence properties can been achieved with `icmod 3`.
In this case, `fcfact` mainly determines the height of the model core charge while
`rcfact` mainly determines the width of the model core charge."""),
        nbf.new_code_cell("fig = plotter.plot_densities(show=False)"),

        nbf.new_markdown_cell("## Local potential and $l$-dependent potentials"),
        nbf.new_code_cell("fig = plotter.plot_potentials(show=False)"),

        #nbf.new_markdown_cell("## 1-st order derivative of $v_l$ and $v_{loc}$ computed via finite differences:"),
        #nbf.new_code_cell("""fig = plotter.plot_der_potentials(order=1, show=False)"""),
        #nbf.new_markdown_cell("## 2-nd order derivative of $v_l$ and $v_{loc}$ computed via finite differences:"),
        #nbf.new_code_cell("""fig = plotter.plot_der_potentials(order=2, show=False)"""),

        nbf.new_markdown_cell("## Model core charge and form factors computed by ABINIT"),
        nbf.new_code_cell("""\
with pseudo.open_pspsfile() as psps:
    fform_fig = psps.plot(show=False);
fform_fig"""),

        nbf.new_markdown_cell("""## Ghosts Test

Self-consistent band structure calculation performed on a regular mesh.
The algorithm to detect ghosts is just an indication usually on the side of false positives.
Zoom in on the band plot to see if an actual ghost is there."""),

        nbf.new_code_cell("fig = report.plot_ebands(with_soc=False, show=False); fig"),

        nbf.new_markdown_cell("""## Convergence of the total energy wrt ecut
The energies are obtained from the deltafactor calculations performed at the Wien2K equilibrium volume"""),
        nbf.new_code_cell("""\
fig = report.plot_etotal_vs_ecut(show=False)"""),

        nbf.new_code_cell("fig = report.plot_etotal_vs_ecut(inv_ecut=True, show=False)"),

        nbf.new_markdown_cell("""## Convergence of the Delta-Gauge results

The Delta-gauge compares the Equation Of State (EOS) of the elemental solid of the element calculated using the pseudo
potential to reference curves calculated using an all electron method. The Delta-gauge was introduced by K. Lejaeghere,
V. Van Speybroeck, G. Van Oost, and&S. Cottenier in [Critical Reviews in Solid State and Materials Sciences 39, 1](http://www.tandfonline.com/doi/abs/10.1080/10408436.2013.772503)

A comparison using the Delta-gauge between many codes and many pseudo tables can be found at the
[center for molecular modeling](http://molmod.ugent.be/deltacodesdft) and in
[Science 351, 1394-1395](http://science.sciencemag.org/content/351/6280/aad3000.full?ijkey=teUZMpwU49vhY&keytype=ref&siteid=sci)
"""),
        nbf.new_code_cell("""fig = report.plot_deltafactor_convergence(xc=pseudo.xc, what=("dfact_meV", "dfactprime_meV"), show=False)"""),
    ])

    # Add validation widget.
    if with_validation:
        nb.cells.extend([
            nbf.new_markdown_cell("## PseudoDojo validation"),
            nbf.new_code_cell("report.ipw_validate()"),
       ])

    nb.cells.extend([
        nbf.new_markdown_cell("## Convergence of $\Delta v_0$, $\Delta b_0$, and $\Delta b_1$ (deltafactor tests)"),
        nbf.new_code_cell("""\
# Absolute difference with respect to Wien2k results.
fig = report.plot_deltafactor_convergence(xc=pseudo.xc, what=("-dfact_meV", "-dfactprime_meV"), show=False)"""),

        nbf.new_markdown_cell("## Delta-gauge EOS for the different cutoff energies"),
        nbf.new_code_cell("fig = report.plot_deltafactor_eos(show=False)"),

        nbf.new_markdown_cell("""## Convergence of the GBRV lattice parameters

The GBRV tests compare the lattice parameter of a FCC and BCC lattice of the element to all electron reference
data. The test was introduced by Kevin F. Garrity, Joseph W. Bennett, Karin M. Rabe, and David Vanderbilt in
developing th GBRV pseudo potential table. More information can be found in [Computational Materials Science 81,
446-452.](http://www.sciencedirect.com/science/article/pii/S0927025613005077)
"""),

        nbf.new_code_cell("fig = report.plot_gbrv_convergence(show=False)"),

        nbf.new_markdown_cell("""## Convergence of the phonon frequencies at $\Gamma$
The calculation is performed with the Wien2k relaxed parameters obtained from the deltafactor CIF files.
"""),
        nbf.new_code_cell("fig = report.plot_phonon_convergence(show=False)"),
    ])

    if with_eos:
        # Add EOS plots
        nb.cells.extend([
            nbf.new_markdown_cell("## GBRV EOS for the FCC structure"),
            nbf.new_code_cell("""fig = report.plot_gbrv_eos(struct_type="fcc", show=False)"""),
            nbf.new_markdown_cell("## GBRV EOS for the BCC structure"),
            nbf.new_code_cell("""fig = report.plot_gbrv_eos(struct_type="bcc", show=False)"""),
        ])

    if tmpfile is None:
        root, ext = os.path.splitext(pseudopath)
        nbpath = root + '.ipynb'
    else:
        import tempfile
        prefix = os.path.basename(pseudopath)
        _, nbpath = tempfile.mkstemp(prefix=prefix, suffix='.ipynb', text=True)

    with io.open(nbpath, 'wt', encoding="utf8") as f:
        nbformat.write(nb, f)

    return nbpath


def make_open_notebook(pseudopath, with_validation=False, with_eos=True):
    """
    Generate a jupyter notebook from the pseudopotential path and
    open it in the browser. Return system exit code.

    Raise:
        RuntimeError if jupyther is not in $PATH
    """
    path = write_notebook(pseudopath, with_validation=with_validation,
                          with_eos=with_eos, tmpfile=True)

    if which("jupyter") is None:
        raise RuntimeError("Cannot find jupyter in PATH. Install it with `pip install`")

    return os.system("jupyter notebook %s" % path)
