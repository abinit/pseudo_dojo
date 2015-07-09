# coding: utf-8
"""Classes and functions for post-processing the results produced by ONCVPSP."""
from __future__ import division, print_function, unicode_literals

import os
import abc
import json
import numpy as np

from collections import namedtuple, OrderedDict
from monty.functools import lazy_property
from monty.collections import AttrDict
from pymatgen.util.plotting_utils import add_fig_kwargs, get_ax_fig_plt
from pseudo_dojo.core import NlState, RadialFunction, RadialWaveFunction
from abipy.tools.derivatives import finite_diff

import logging
logger = logging.getLogger(__name__)


_l2char = {
    "0": "s",
    "1": "p",
    "2": "d",
    "3": "f",
    "4": "g",
    "5": "h",
    "6": "i",
}

def is_integer(s):
    try:
        c = float(s)
        return int(c) == c
    except (ValueError, TypeError):
        return False


def decorate_ax(ax, xlabel, ylabel, title, lines, legends):
    """Decorate a `matplotlib` Axis adding xlabel, ylabel, title, grid and legend"""
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True)
    ax.legend(lines, legends, loc="best", shadow=True)



class PseudoGenDataPlotter(object):
    """
    Plots the results produced by a pseudopotential generator.
    """
    # TODO: Add fully-relativistic case.
    # List of results supported by the plotter (initialized via __init__)
    all_keys = [
        "radial_wfs",
        "projectors",
        "densities",
        "potentials",
        "atan_logders",
        "ene_vs_ecut",
    ]

    # matplotlib options.
    linewidth, markersize = 2, 2

    linestyle_aeps = dict(ae="solid", ps="dashed")
    markers_aeps = dict(ae=".", ps="o")
    color_l = {-1: "black", 0: "red", 1: "blue", 2: "green", 3: "orange"}

    def __init__(self, **kwargs):
        """Store kwargs in self if k is in self.all_keys."""
        import matplotlib.pyplot as _mplt
        self._mplt = _mplt

        for k in self.all_keys:
            setattr(self, k, kwargs.pop(k, {}))

        if kwargs:
            raise ValueError("Unknown keys: %s" % list(kwargs.keys()))

    def keys(self):
        """Iterator with the keys stored in self."""
        return (k for k in self.all_keys if getattr(self, k))

    def plot_key(self, key, ax=None, **kwargs):
        """Plot a singol quantity specified by key."""
        ax, fig, plt = get_ax_fig_plt(ax)

        # key --> self.plot_key()
        getattr(self, "plot_" + key)(ax=ax, **kwargs)
        self._mplt.show()

    def _wf_pltopts(self, l, aeps):
        """Plot options for wavefunctions."""
        return dict(
            color=self.color_l[l], linestyle=self.linestyle_aeps[aeps], #marker=self.markers_aeps[aeps],
            linewidth=self.linewidth, markersize=self.markersize)

    @add_fig_kwargs
    def plot_atan_logders(self, ax=None, **kwargs):
        """Plot arctan of logder on axis ax."""
        ae, ps = self.atan_logders.ae, self.atan_logders.ps

        ax, fig, plt = get_ax_fig_plt(ax)

        lines, legends = [], []
        for l, ae_alog in ae.items():
            ps_alog = ps[l]

            # Add padd to avoid overlapping curves.
            pad = (l+1) * 1.0

            ae_line, = ax.plot(ae_alog.energies, ae_alog.values + pad, **self._wf_pltopts(l, "ae"))
            ps_line, = ax.plot(ps_alog.energies, ps_alog.values + pad, **self._wf_pltopts(l, "ps"))

            lines.extend([ae_line, ps_line])
            legends.extend(["AE l=%s" % str(l), "PS l=%s" % str(l)])

        decorate_ax(ax, xlabel="Energy [Ha]", ylabel="ATAN(LogDer)", title="ATAN(Log Derivative)", 
                    lines=lines, legends=legends)

        return fig

    @add_fig_kwargs
    def plot_radial_wfs(self, ax=None, **kwargs):
        """
        Plot ae and ps radial wavefunctions on axis ax.

        lselect: List to select l channels
        """
        ax, fig, plt = get_ax_fig_plt(ax)

        ae_wfs, ps_wfs = self.radial_wfs.ae, self.radial_wfs.ps
        lselect = kwargs.get("lselect", [])

        lines, legends = [], []
        for nl, ae_wf in ae_wfs.items():
            ps_wf, l = ps_wfs[nl], nl.l
            if l in lselect: continue

            ae_line, = ax.plot(ae_wf.rmesh, ae_wf.values, **self._wf_pltopts(l, "ae"))
            ps_line, = ax.plot(ps_wf.rmesh, ps_wf.values, **self._wf_pltopts(l, "ps"))

            lines.extend([ae_line, ps_line])
            legends.extend(["AE l=%s" % str(l), "PS l=%s" % str(l)])

        decorate_ax(ax, xlabel="r [Bohr]", ylabel="$\phi(r)$", title="Wave Functions", 
                    lines=lines, legends=legends)

        return fig

    @add_fig_kwargs
    def plot_projectors(self, ax=None, **kwargs):
        """
        Plot oncvpsp projectors on axis ax.

        lselect: List to select l channels
        """
        ax, fig, plt = get_ax_fig_plt(ax)

        lselect = kwargs.get("lselect", [])

        lines, legends = [], []
        for nl, proj in self.projectors.items():
            if nl.l in lselect: continue
            line, = ax.plot(proj.rmesh, proj.values,
                            linewidth=self.linewidth, markersize=self.markersize)

            lines.append(line)
            legends.append("Proj %s" % str(nl))

        decorate_ax(ax, xlabel="r [Bohr]", ylabel="$p(r)$", title="Projector Wave Functions", 
                    lines=lines, legends=legends)
        return fig

    @add_fig_kwargs
    def plot_densities(self, ax=None, timesr2=False, **kwargs):
        """Plot ae, ps and model densities on axis ax."""
        ax, fig, plt = get_ax_fig_plt(ax)

        lines, legends = [], []
        for name, rho in self.densities.items():
            d = rho.values if not timesr2 else  rho.values * rho.rmesh ** 2 
            line, = ax.plot(rho.rmesh, d,
                            linewidth=self.linewidth, markersize=self.markersize)

            lines.append(line)
            legends.append(name)

        ylabel = "$n(r)$" if not timesr2 else "$r^2 n(r)$"
        decorate_ax(ax, xlabel="r [Bohr]", ylabel=ylabel, title="Charge densities", 
                    lines=lines, legends=legends)

        return fig

    @add_fig_kwargs
    def plot_der_densities(self, ax=None, order=1, **kwargs):
        """
        Plot the derivatives of the densitiers on axis ax.
        Used to analyze possible derivative discontinuities
        """
        ax, fig, plt = get_ax_fig_plt(ax)

        from scipy.interpolate import UnivariateSpline

        lines, legends = [], []
        for name, rho in self.densities.items():
            # Need linear mesh for finite_difference --> Spline input densities on lin_rmesh
            lin_rmesh, h = np.linspace(rho.rmesh[0], rho.rmesh[-1], num=len(rho.rmesh) * 4, retstep=True)
            spline = UnivariateSpline(rho.rmesh, rho.values, s=0)
            lin_values = spline(lin_rmesh)
            vder = finite_diff(lin_values, h, order=order, acc=4)
            line, = ax.plot(lin_rmesh, vder) #, **self._wf_pltopts(l, "ae"))
            lines.append(line)
                                                                                             
            legends.append("%s-order derivative of %s" % (order, name))
                                                                                             
        decorate_ax(ax, xlabel="r [Bohr]", ylabel="$D^%s \n(r)$" % order, title="Derivative of the charge densities", 
                    lines=lines, legends=legends)
        return fig

    @add_fig_kwargs
    def plot_potentials(self, ax=None, **kwargs):
        """Plot vl and vloc potentials on axis ax"""
        ax, fig, plt = get_ax_fig_plt(ax)

        lines, legends = [], []
        for l, pot in self.potentials.items():
            line, = ax.plot(pot.rmesh, pot.values, **self._wf_pltopts(l, "ae"))
            lines.append(line)

            if l == -1:
                legends.append("Vloc")
            else:
                legends.append("PS l=%s" % str(l))

        decorate_ax(ax, xlabel="r [Bohr]", ylabel="$v_l(r)$", title="Ion Pseudopotentials", 
                    lines=lines, legends=legends)
        return fig

    @add_fig_kwargs
    def plot_der_potentials(self, ax=None, order=1, **kwargs):
        """
        Plot the derivatives of vl and vloc potentials on axis ax.
        Used to analyze the derivative discontinuity introduced by the RRKJ method at rc.
        """
        ax, fig, plt = get_ax_fig_plt(ax)
        from abipy.tools.derivatives import finite_diff
        from scipy.interpolate import UnivariateSpline
        lines, legends = [], []
        for l, pot in self.potentials.items():
            # Need linear mesh for finite_difference --> Spline input potentials on lin_rmesh
            lin_rmesh, h = np.linspace(pot.rmesh[0], pot.rmesh[-1], num=len(pot.rmesh) * 4, retstep=True)
            spline = UnivariateSpline(pot.rmesh, pot.values, s=0)
            lin_values = spline(lin_rmesh)
            vder = finite_diff(lin_values, h, order=order, acc=4)
            line, = ax.plot(lin_rmesh, vder, **self._wf_pltopts(l, "ae"))
            lines.append(line)
                                                                                             
            if l == -1:
                legends.append("%s-order derivative Vloc" % order)
            else:
                legends.append("$s-order derivative PS l=%s" % str(l))
                                                                                             
        decorate_ax(ax, xlabel="r [Bohr]", ylabel="$D^%s \phi(r)$" % order, title="Derivative of the ion Pseudopotentials", 
                    lines=lines, legends=legends)
        return fig

    @add_fig_kwargs
    def plot_ene_vs_ecut(self, ax=None, **kwargs):
        """Plot the converge of ene wrt ecut on axis ax."""
        ax, fig, plt = get_ax_fig_plt(ax)
        lines, legends = [], []
        for l, data in self.ene_vs_ecut.items():
            line, = ax.plot(data.energies, data.values, **self._wf_pltopts(l, "ae"))

            lines.append(line)
            legends.append("Conv l=%s" % str(l))

        decorate_ax(ax, xlabel="Ecut [Ha]", ylabel="$\Delta E$", title="Energy error per electron [Ha]", 
                    lines=lines, legends=legends)

        ax.set_yscale("log")
        return fig

    @add_fig_kwargs
    def plot_atanlogder_econv(self, **kwargs):
        """Plot atan(logder) and ecut converge on the same figure. Returns matplotlib Figure"""
        fig, ax_list = self._mplt.subplots(nrows=2, ncols=1, sharex=False, squeeze=False)
        ax_list = ax_list.ravel()

        self.plot_atan_logders(ax=ax_list[0], show=False)
        self.plot_ene_vs_ecut(ax=ax_list[1], show=False)

        return fig

    @add_fig_kwargs
    def plot_dens_and_pots(self, **kwargs):
        """Plot densities and potentials on the same figure. Returns matplotlib Figure"""
        fig, ax_list = self._mplt.subplots(nrows=2, ncols=1, sharex=False, squeeze=False)
        ax_list = ax_list.ravel()

        self.plot_densities(ax=ax_list[0], show=False)
        self.plot_potentials(ax=ax_list[1], show=False)

        return fig

    @add_fig_kwargs
    def plot_waves_and_projs(self, **kwargs):
        """Plot ae-ps wavefunctions and projectors on the same figure. Returns matplotlib Figure"""
        lmax = max(nl.l for nl in self.radial_wfs.ae.keys())
        fig, ax_list = self._mplt.subplots(nrows=lmax+1, ncols=2, sharex=True, squeeze=False)

        for l in range(lmax+1):
            ax_idx = lmax - l
            self.plot_radial_wfs(ax=ax_list[ax_idx][0], lselect=[l], show=False)
            self.plot_projectors(ax=ax_list[ax_idx][1], lselect=[l], show=False)

        return fig

    @add_fig_kwargs
    def plot_den_formfact(self, ecut=60, ax=None, **kwargs):
        """
        Plot the density form factor as function of ecut (Ha units). Return matplotlib Figure.
        """
        ax, fig, plt = get_ax_fig_plt(ax)

        lines, legends = [], []
        for name, rho in self.densities.items():
            if name == "rhoC": continue
            form = rho.get_intr2j0(ecut=ecut) / (4 * np.pi)
            line, = ax.plot(form.mesh, form.values, linewidth=self.linewidth, markersize=self.markersize)
            lines.append(line); legends.append(name)

            intg = rho.r2f_integral()[-1]
            print("r2 f integral: ", intg)
            print("form_factor(0): ", name, form.values[0]) 

        # Plot vloc(q)
        #for l, pot in self.potentials.items():
        #    if l != -1: continue
        #    form = pot.get_intr2j0(ecut=ecut)
        #    mask = np.where(np.abs(form.values) > 20); form.values[mask] = 20
        #    line, = ax.plot(form.mesh, form.values, linewidth=self.linewidth, markersize=self.markersize)
        #    lines.append(line); legends.append("Vloc(q)")

        decorate_ax(ax, xlabel="Ecut [Ha]", ylabel="$n(q)$", title="Form factor, l=0 ", lines=lines, legends=legends)
        return fig


class MultiPseudoGenDataPlotter(object):
    """
    Class for plotting data produced by multiple pp generators on separated plots.

    Usage example:

    .. code-block:: python

        plotter = MultiPseudoGenPlotter()
        plotter.add_psgen("bar.nc", label="bar bands")
        plotter.add_plotter("foo.nc", label="foo bands")
        plotter.plot()
    """
    _LINE_COLORS = ["b", "r"]
    _LINE_STYLES = ["-", ":", "--", "-."]
    _LINE_WIDTHS = [2]

    def __init__(self):
        self._plotters_odict = OrderedDict()

    def __len__(self):
        return len(self._plotters_odict)

    @property
    def plotters(self):
        """"List of registered `Plotters`."""
        return list(self._plotters_odict.values())

    @property
    def labels(self):
        """List of labels."""
        return list(self._plotters_odict.keys())

    def keys(self):
        """List of strings with the quantities that can be plotted."""
        keys_set = set()
        for plotter in self.plotters:
            keys_set.update(plotter.keys())
        keys_set = list(keys_set)

        return list(sorted(keys_set))

    def iter_lineopt(self):
        """Generates style options for lines."""
        import itertools
        for o in itertools.product(self._LINE_WIDTHS,  self._LINE_STYLES, self._LINE_COLORS):
            yield {"linewidth": o[0], "linestyle": o[1], "color": o[2]}

    def add_psgen(self, label, psgen):
        """Add a plotter of class plotter_class from a `PseudoGenerator` instance."""
        oparser = psgen.parse_output()
        self.add_plotter(label, oparser.make_plotter())

    def add_plotter(self, label, plotter):
        """
        Adds a plotter.

        Args:
            label: label for the plotter. Must be unique.
            plotter: :class:`PseudoGenDataPlotter` object.
        """
        if label in self.labels:
            raise ValueError("label %s is already in %s" % (label, self.labels))

        self._plotters_odict[label] = plotter

    @add_fig_kwargs
    def plot_key(self, key, **kwargs):
        """
        Plot the band structure and the DOS.

        Args:
            klabels: dictionary whose keys are tuple with the reduced
                coordinates of the k-points. The values are the labels.
                e.g. klabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0): "L"}.

        ==============  ==============================================================
        kwargs          Meaning
        ==============  ==============================================================
        sharex          True if subplots should share the x axis
        ==============  ==============================================================

        Returns:
            matplotlib figure.
        """
        import matplotlib.pyplot as plt

        # Build grid of plots.
        fig, ax_list = plt.subplots(nrows=len(self), ncols=1, sharex=kwargs.pop("sharex", True), squeeze=False)
        ax_list = ax_list.ravel()

        for ax in ax_list:
            ax.grid(True)

        # Plot key on the each axis in ax_list.
        lines, legends = [], []
        #my_kwargs, opts_label = kwargs.copy(), {}
        i = -1
        for (label, plotter), lineopt in zip(self._plotters_odict.items(), self.iter_lineopt()):
            i += 1
            plotter.plot_key(key, ax=ax_list[i])

        return fig


class PseudoGenResults(AttrDict):
    _KEYS = [
        "max_ecut",
        "max_atan_logder_l1err",
    ]

    def __init__(self, *args, **kwargs):
        super(PseudoGenResults, self).__init__(*args, **kwargs)
        for k in self._KEYS:
            if k not in self:
                self[k] = None


class AtanLogDer(namedtuple("AtanLogDer", "l, energies, values")):
    @property
    def to_dict(self):
        return dict(
            l=self.l,
            energies=list(self.energies),
            values=list(self.values))


class PseudoGenOutputParserError(Exception):
    """Exceptions raised by OuptputParser."""


class PseudoGenOutputParser(object):
    """
    Abstract class defining the interface that must be provided
    by the parsers used to extract results from the output file of
    a pseudopotential generator a.k.a. ppgen

    Attributes:

        errors: List of strings with errors reported by the pp generator
        warnings: List of strings with the warnings reported by the pp generator.
        results:
    """
    __metaclass__ = abc.ABCMeta

    Error = PseudoGenOutputParserError

    def __init__(self, filepath):
        self.filepath = os.path.abspath(filepath)
        self.run_completed = False
        self._errors = []
        self._warnings = []
        self._results = None

    @property
    def errors(self):
        """List of strings with possible errors reported by the generator at run-time."""
        return self._errors

    @property
    def warnings(self):
        """List of strings with possible errors reported by the generator at run-time."""
        return self._warnings

    @property
    def results(self):
        return self._results

    @abc.abstractmethod
    def get_results(self):
        """
        Return the most important results of the run in a dictionary.
        Set self.results attribute
        """

    @abc.abstractmethod
    def get_input_str(self):
        """Returns a string with the input file."""
    
    @abc.abstractmethod
    def get_pseudo_str(self):
        """Returns a string with the pseudopotential file."""


class OncvOutputParser(PseudoGenOutputParser):
    """
    Object to read and extract data from the output file of oncvpsp.

    Attributes:
        atsym
        Z
        nc
        nv
        iexc
        psfile

    Example:
        parser = OncvOutputParser(filename)

        # To access data:
        p.radial_wavefunctions

        # To plot data with matplotlib.
        p = parser.make_plotter()
        p.plot_slideshow()
    """
    # TODO Add fully-relativistic case.

    # Used to store ae and pp quantities (e.g wavefunctions) in a single object.
    AePsNamedTuple = namedtuple("AePsNamedTuple", "ae, ps")

    # Object returned by self._grep
    GrepResults = namedtuple("GrepResults", "data, start, stop")

    # Object storing the final results.
    Results = PseudoGenResults

    # Class used to instanciate the plotter.
    Plotter = PseudoGenDataPlotter

    def scan(self):
        """
        Scan the output, set and set `run_completed` attribute.

        Raises:
            self.Error if invalid file.
        """
        if not os.path.exists(self.filepath):
            raise self.Error("File %s does not exist" % self.filepath)

        # Read data and store it in lines
        self.lines = []
        with open(self.filepath) as fh:
            for i, line in enumerate(fh):
                line = line.strip()
                self.lines.append(line)

                if line.startswith("DATA FOR PLOTTING"):
                    self.run_completed = True

                if "ERROR" in line:
                    # Example:
                    # test_data: must have fcfact>0.0 for icmod= 1
                    # ERROR: test_data found   1 errors; stopping
                    self._errors.append("\n".join(self.lines[i-1:i+1]))

                if "WARNING" in line:
                    self._warnings.append("\n".join(self.lines[i:i+2]))

        #if self.errors:
        #    return 1

        #scalar-relativistic version 2.1.1, 03/26/2014
        #scalar-relativistic version 3.0.0 10/10/2014
        #toks, self.gendate = self.lines[1].split(",")
        #toks = toks.split()
        toks = self.lines[1].replace(",", " ").split()
        self.gendate = toks.pop(-1)
        self.calc_type, self.version = toks[0], toks[-1]

        if self.calc_type not in ["scalar-relativistic", "non-relativistic"]:
            print("will raise %s because found %s" % (self.Error, self.calc_type))
            raise self.Error("Fully relativistic case is not supported")

        # Read configuration (not very robust because we assume the user didn't change the template but oh well)
        header = "# atsym  z    nc    nv    iexc   psfile"
        for i, line in enumerate(self.lines):
            if line.startswith("# atsym"):
                values = self.lines[i+1].split()
                keys = header[1:].split()
                assert len(keys) == len(values)
                # Store them in self.
                for k, v in zip(keys, values):
                    setattr(self, k, v)
                break

        # Parse ATOM and Rerencence configuration
        # Example::
        """
        #
        #   n    l    f        energy (Ha)
            1    0    2.00    -6.5631993D+01
            2    0    2.00    -5.1265474D+00
            2    1    6.00    -3.5117357D+00
            3    0    2.00    -3.9736459D-01
            3    1    2.00    -1.4998149D-01
        """
        header = "#   n    l    f        energy (Ha)"
        nc, nv = int(self.nc), int(self.nv)
        for i, line in enumerate(self.lines):
            if line.startswith(header):
                beg, core = i + 1, [], 
                for c in range(nc):
                    n, l, f = self.lines[beg+c].split()[:3]
                    if is_integer(f):
                        f = str(int(float(f)))
                    else:
                        f = "%.1f" % f
                    core.append(n + _l2char[l] + "^%s" %f)
                self.core = "$" + " ".join(core) + "$"

                beg, valence = i + nc + 1, [] 
                for v in range(nv):
                    n, l, f = self.lines[beg+v].split()[:3]
                    if is_integer(f):
                        f = str(int(float(f)))
                    else:
                        f = "%.1f" % f

                    valence.append(n + _l2char[l] + "^{%s}" % f)
                self.valence = "$" + " ".join(valence) + "$"
                #print("core", self.core)
                #print("valence",self.valence)
                break
        else:
            raise self.Error("Cannot find #lmax line in output file %s" % self.filepath)

        # Read lmax (not very robust because we assume the user didn't change the template but oh well)
        header = "# lmax"
        for i, line in enumerate(self.lines):
            if line.startswith(header):
                self.lmax = int(self.lines[i+1])
                break
        else:
            raise self.Error("Cannot find #lmax line in output file %s" % self.filepath)
        #print("lmax", self.lmax)

        # Compute the minimun rc(l)
        header = "#   l,   rc,     ep,   ncon, nbas, qcut"
        for i, line in enumerate(self.lines):
            if line.startswith(header):
                beg = i + 1
                next, rcs = 0, []
                while True:
                    l = self.lines[beg + next]
                    if l.startswith("#"): break
                    token = l.split()[1]
                    #print("token", token)
                    rcs.append(float(token))
                    #print(l)
                    next += 1

                self.rc_min, self.rc_max = min(rcs), max(rcs)

    def __str__(self):
        """String representation."""
        lines = []
        app = lines.append
        app("%s, oncvpsp version: %s, date: %s" % (self.calc_type, self.version, self.gendate))
        app("oncvpsp calculation: %s: " % self.calc_type)
        app("completed: %s" % self.run_completed)

        return "\n".join(lines)

    @property
    def fully_relativistic(self):
        """True if fully-relativistic calculation."""
        return self.calc_type == "fully-relativistic"

    @lazy_property
    def potentials(self):
        """Radial functions with the non-local and local potentials."""
        #radii, charge, pseudopotentials (ll=0, 1, lmax)
        #!p   0.0099448   4.7237412  -7.4449470 -14.6551019
        vl_data = self._grep("!p").data
        lmax = len(vl_data[0]) - 3
        assert lmax == self.lmax

        # From 0 up to lmax
        ionpots_l = {}
        for l in range(lmax+1):
            ionpots_l[l] = RadialFunction("Ion Pseudopotential, l=%d" % l, vl_data[:, 0], vl_data[:, 2+l])

        # Local part is stored with l == -1 if lloc=4, not present if lloc=l
        vloc = self._grep("!L").data
        if vloc is not None:
            ionpots_l[-1] = RadialFunction("Local part, l=%d" % -1, vloc[:, 0], vloc[:, 1])

        return ionpots_l

    @lazy_property
    def densities(self):
        """Dictionary with charge densities on the radial mesh."""
        # radii, charge, core charge, model core charge
        # !r   0.0100642   4.7238866  53.4149287   0.0000000
        rho_data = self._grep("!r").data

        return dict(
            rhoV=RadialFunction("Valence charge", rho_data[:, 0], rho_data[:, 1]),
            rhoC=RadialFunction("Core charge", rho_data[:, 0], rho_data[:, 2]),
            rhoM=RadialFunction("Model charge", rho_data[:, 0], rho_data[:, 3]))

    @lazy_property
    def radial_wfs(self):
        """Read the radial wavefunctions."""
        #n= 1,  l= 0, all-electron wave function, pseudo w-f
        #
        #&     0    0.009945   -0.092997    0.015273
        ae_waves, ps_waves = OrderedDict(), OrderedDict()

        beg = 0
        while True:
            g = self._grep("&", beg=beg)
            if g.data is None:
                break
            beg = g.stop + 1

            header = self.lines[g.start-2]
            n, l = header.split(",")[0:2]
            n = int(n.split("=")[1])
            l = int(l.split("=")[1])
            nl = NlState(n=n, l=l)
            logger.info("Got state: %s" % str(nl))

            rmesh = g.data[:, 1]
            ae_wf = g.data[:, 2]
            ps_wf = g.data[:, 3]

            ae_waves[nl] = RadialWaveFunction(nl, str(nl), rmesh, ae_wf)
            ps_waves[nl] = RadialWaveFunction(nl, str(nl), rmesh, ps_wf)

        return self.AePsNamedTuple(ae=ae_waves, ps=ps_waves)

    @lazy_property
    def projectors(self):
        """Read the projector wave functions."""
        #n= 1 2  l= 0, projecctor pseudo wave functions, well or 2nd valence
        #
        #@     0    0.009945    0.015274   -0.009284
        projectors_nl = OrderedDict()
        beg = 0
        while True:
            g = self._grep("@", beg=beg)
            if g.data is None:
                break
            beg = g.stop + 1

            rmesh = g.data[:, 1]
            l = int(g.data[0, 0])

            for n in range(len(g.data[0]) - 2):
                nl = NlState(n=n+1, l=l)
                logger.info("Got projector with: %s" % str(nl))
                projectors_nl[nl] = RadialWaveFunction(nl, str(nl), rmesh, g.data[:, n+2])

        return projectors_nl

    @lazy_property
    def atan_logders(self):
        """Atan of the log derivatives for different l-values."""
        #log derivativve data for plotting, l= 0
        #atan(r * ((d psi(r)/dr)/psi(r))), r=  1.60
        #l, energy, all-electron, pseudopotential
        #
        #!      0    2.000000    0.706765    0.703758
        ae_atan_logder_l, ps_atan_logder_l = OrderedDict(), OrderedDict()

        for l in range(self.lmax+1):
            data = self._grep(tag="!      %d" % l).data
            assert l == int(data[0, 0])
            ae_atan_logder_l[l] = AtanLogDer(l=l, energies=data[:, 1], values=data[:, 2])
            ps_atan_logder_l[l] = AtanLogDer(l=l, energies=data[:, 1], values=data[:, 3])

        return self.AePsNamedTuple(ae=ae_atan_logder_l, ps=ps_atan_logder_l)

    @lazy_property
    def ene_vs_ecut(self):
        """Convergence of energy versus ecut for different l-values."""
        #convergence profiles, (ll=0,lmax)
        #!C     0    5.019345    0.010000
        #...
        #!C     1   19.469226    0.010000
        class ConvData(namedtuple("ConvData", "l energies values")):
            @property
            def to_dict(self):
                return dict(
                    l=self.l,
                    energies=list(self.energies),
                    values=list(self.values))

        conv_l = OrderedDict()

        for l in range(self.lmax+1):
            data = self._grep(tag="!C     %d" % l).data
            conv_l[l] = ConvData(l=l, energies=data[:, 1], values=data[:, 2])

        return conv_l

    @lazy_property
    def hints(self):
        """Hints for the cutoff energy."""
        # Extract the hints
        hints = 3 * [-np.inf]
        ene_vs_ecut = self.ene_vs_ecut
        for i in range(3):
            for l in range(self.lmax+1):
                hints[i] = max(hints[i], ene_vs_ecut[l].energies[-i-1])
        hints.reverse()

        # print("hints:", hints)
        # Truncate to the nearest int
        hints = [np.rint(h) for h in hints]

        hints = dict(
            low={"ecut": hints[0], "pawecutdg": hints[0]},
            normal={"ecut": hints[1], "pawecutdg": hints[1]},
            high={"ecut": hints[2], "pawecutdg": hints[2]})

        return hints

    def get_results(self):
        """"
        Return the most important results reported by the pp generator.
        Set the value of self.results
        """
        #if not self.run_completed:
        #    self.Results(info="Run is not completed")

        # Get the ecut needed to converge within ... TODO
        max_ecut = 0.0
        for l in range(self.lmax+1):
            max_ecut = max(max_ecut, self.ene_vs_ecut[l].energies[-1])

        # Compute the l1 error in atag(logder)
        from scipy.integrate import cumtrapz
        max_l1err = 0.0
        for l in range(self.lmax+1):
            f1, f2 = self.atan_logders.ae[l], self.atan_logders.ps[l]

            adiff = np.abs(f1.values - f2.values)
            integ = cumtrapz(adiff, x=f1.energies) / (f1.energies[-1] - f1.energies[0])
            max_l1err = max(max_l1err, integ[-1])

        # Read Hermiticity error and compute the max value of PSP excitation error=
        # Hermiticity error    4.8392D-05
        # PSP excitation error=  1.56D-10
        herm_tag, pspexc_tag = "Hermiticity error", "PSP excitation error="
        herm_err, max_psexc_abserr = None, -np.inf

        for line in self.lines:
            i = line.find(herm_tag)
            if i != -1:
                herm_err = float(line.split()[-1].replace("D", "E"))

            i = line.find(pspexc_tag)
            if i != -1:
                max_psexc_abserr = max(max_psexc_abserr, abs(float(line.split()[-1].replace("D", "E"))))

        self._results = self.Results(
            max_ecut=max_ecut, max_atan_logder_l1err=max_l1err,
            herm_err=herm_err, max_psexc_abserr=max_psexc_abserr)

        return self._results

    def find_string(self, s):
        """
        Returns the index of the first line containing string s.
        Raise self.Error if s cannot be found.
        """
        for i, line in enumerate(self.lines):
            if s in line:
                return i
        else:
            raise self.Error("Cannot find %s in lines" % s)

    def get_input_str(self):
        """String with the input file."""
        i = self.find_string("Reference configufation results")
        return "\n".join(self.lines[:i])

    def get_pseudo_str(self, devel=False):
        """String with the pseudopotential data."""
        # devel is for tuning the pseudo, only two cutoffs

        # Extract the pseudo in Abinit format.
        try:
            i = self.find_string('Begin PSPCODE8')
        except:
            i = self.find_string('Begin PSP_UPF')

        ps_data = "\n".join(self.lines[i+1:])

        # Append the input to ps_data (note XML markers)
        ps_data += "\n\n<INPUT>\n" + self.get_input_str() + "</INPUT>\n"

        # Add the initial DOJO_REPORT with the hints and the initial list of ecut values.
        estart = self.hints["high"]["ecut"]
        dense_right = np.linspace(estart - 10, estart + 10, num=11)

        d = {
            "version": "1.0",
            "pseudo_type": "norm-conserving",
            "ppgen_hints": self.hints, 
            "ecuts": list(dense_right) + [dense_right[-1] + 8, dense_right[-1] + 10,],
            "symbol": self.atsym,
        }

        if devel:
            # devellopment run: few, relatively high ecut calculations
            d["ecuts"] = [estart - 10, estart + 10]

        ps_data += "\n<DOJO_REPORT>\n" + json.dumps(d, indent=4) + "\n</DOJO_REPORT>\n"

        return ps_data

    def make_plotter(self):
        """Builds an instance of :class:`PseudoGenDataPlotter`."""
        kwargs = {k: getattr(self, k) for k in self.Plotter.all_keys}
        return self.Plotter(**kwargs)

    @property
    def to_dict(self):
        """
        Returns a dictionary with the radial functions and the other
        important results produced by ONCVPSP in JSON format.
        """
        # Dimensions and basic info.
        jdict = dict(
            lmax=self.lmax,
            ncore=self.nc,
            nvalence=self.nv,
            calc_type=self.calc_type)

        # List of radial wavefunctions (AE and PS)
        jdict["radial_wfs"] = d = {}
        d["ae"] = [wave.to_dict for wave in self.radial_wfs.ae.values()]
        d["ps"] = [wave.to_dict for wave in self.radial_wfs.ps.values()]

        # List of projectors
        jdict["projectors"] = [proj.to_dict for proj in self.projectors.values()]

        # Charge densities
        jdict["densities"] = dict(
            rhoV=self.densities["rhoV"].to_dict,
            rhoC=self.densities["rhoC"].to_dict,
            rhoM=self.densities["rhoM"].to_dict)

        # Logders (AE and PS)
        jdict["atan_logders"] = d = {}
        d["ae"] = [f.to_dict for f in self.atan_logders.ae.values()]
        d["ps"] = [f.to_dict for f in self.atan_logders.ps.values()]

        # Convergence of the different l-channels as function of ecut.
        jdict["ene_vs_ecut"] = [f.to_dict for f in self.ene_vs_ecut.values()]

        return jdict

    def _grep(self, tag, beg=0):
        """
        Finds the first field in the file with the specified tag.
        beg gives the initial position in the file.

        Returns:
            :class:`GrepResult` object
        """
        data, stop, intag = [], None, -1

        if beg >= len(self.lines):
            raise ValueError()

        for i, l in enumerate(self.lines[beg:]):
            l = l.lstrip()
            if l.startswith(tag):
                if intag == -1:
                    intag = beg + i
                data.append([float(c) for c in l.split()[1:]])
            else:
                # Exit because we know there's only one section starting with 'tag'
                if intag != -1:
                    stop = beg + i
                    break

        if not data:
            return self.GrepResults(data=None, start=intag, stop=stop)
        else:
            return self.GrepResults(data=np.array(data), start=intag, stop=stop)
