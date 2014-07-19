#!/usr/bin/env python
"""Classes and functions for post-processing the results produced by ONCVPSP."""
from __future__ import print_function, division

import os
import collections
import numpy as np

from pprint import pprint
from pseudo_dojo.core import NlState, RadialFunction, RadialWaveFunction


class OncvOutputParserError(Exception):
    """Exceptions raised by OncvParser."""


class OncvOuptputParser(object):
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
    Error = OncvOutputParserError

    # Used to store ae and pp quantities (e.g wavefunctions) in a single object.
    AePsNamedTuple = collections.namedtuple("AePsNamedTuple", "ae, ps")

    # object returned by self._grep
    GrepResults = collections.namedtuple("GrepResults", "data, start, stop")

    def __init__(self, filepath):
        """Initialize the object from the oncvpsp output."""
        self.filepath = os.path.abspath(filepath)
        # Read data and store it in lines
        with open(self.filepath) as fh:
            self.lines = fh.readlines()

        # scalar-relativistic version 2.1.1, 03/26/2014
        toks, self.gendate = self.lines[1].split(",")

        toks = toks.split()
        self.calc_type, self.version = toks[0], toks[-1]

        if self.calc_type not in ["scalar-relativistic", "non-relativistic"]:
            raise self.Error("Fully relativistic case is not supported")

        # Read configuration (not very robust because we assume the user didn't change the template but oh well)
        header = "# atsym  z    nc    nv    iexc   psfile"
        for i, line in enumerate(self.lines):
            if line.startswith(header):
                values = self.lines[i+1].split()
                keys = header[1:].split()
                assert len(keys) == len(values)
                # Store them in self.
                for k, v in zip(keys, values):
                    setattr(self, k, v)
                break

        # Read lmax (not very robust because we assume the user didn't change the template but oh well)
        header = "# lmax"
        for i, line in enumerate(self.lines):
            if line.startswith(header):
                self.lmax = int(self.lines[i+1])
                break
        else:
            raise self.Error("Cannot find #lmax line in output file %s" % self.filepath)

    def __str__(self):
        lines = []
        app = lines.append
        app("%s, oncvpsp version: %s, date: %s" % (self.calc_type, self.version, self.gendate))
        return "\n".join(lines)

    @property
    def fully_relativistic(self):
        return self.calc_type == "fully-relativistic"

    @property
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

        # Local part is stored with l == -1
        vloc = self._grep("!L").data
        ionpots_l[-1] = RadialFunction("Local part, l=%d" % -1, vloc[:, 0], vloc[:, 1])

        return ionpots_l

    @property
    def densities(self):
        """Dictionary with charge densities on the radial mesh."""
        # radii, charge, core charge, model core charge
        # !r   0.0100642   4.7238866  53.4149287   0.0000000
        rho_data = self._grep("!r").data

        return dict(
            rhoV=RadialFunction("Valence charge", rho_data[:, 0], rho_data[:, 1]),
            rhoC=RadialFunction("Core charge", rho_data[:, 0], rho_data[:, 2]),
            rhoM=RadialFunction("Model charge", rho_data[:, 0], rho_data[:, 3]))

    @property
    def radial_wfs(self):
        """Read the radial wavefunctions."""
        # TODO: Check this
        ae_waves, ps_waves = {}, {}
        #n= 1,  l= 0, all-electron wave function, pseudo w-f
        #
        #&     0    0.009945   -0.092997    0.015273

        beg = 0
        while True:
            g = self._grep("&", beg=beg)
            #print(g)
            if g.data is None:
                break
            beg = g.stop + 1

            header = self.lines[g.start-2]
            n, l = header.split(",")[0:2]
            n = int(n.split("=")[1])
            l = int(l.split("=")[1])
            nl = NlState(n=n, l=l)
            #print("Got state: %s" % str(nl))

            rmesh = g.data[:, 1]
            ae_wf = g.data[:, 2]
            ps_wf = g.data[:, 3]

            ae_waves[nl] = RadialWaveFunction(nl, str(nl), rmesh, ae_wf)
            ps_waves[nl] = RadialWaveFunction(nl, str(nl), rmesh, ps_wf)

        return self.AePsNamedTuple(ae=ae_waves, ps=ps_waves)

    @property
    def projectors(self):
        """Read the projector wave functions."""
        # TODO: Check this
        #n= 1 2  l= 0, projecctor pseudo wave functions, well or 2nd valence
        #
        #@     0    0.009945    0.015274   -0.009284
        projectors_nl = {}
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
                #print("Got projector with: %s" % str(nl))
                projectors_nl[nl] = RadialWaveFunction(nl, str(nl), rmesh, g.data[:, n+2])

        return projectors_nl

    @property
    def atan_logders(self):
        """Atan of the log derivatives for different l-values."""
        #log derivativve data for plotting, l= 0
        #atan(r * ((d psi(r)/dr)/psi(r))), r=  1.60
        #l, energy, all-electron, pseudopotential
        #
        #!      0    2.000000    0.706765    0.703758
        atan_logder = collections.namedtuple("AtanLogDer", "energies values")
        ae_atan_logder_l, ps_atan_logder_l = {}, {}

        for l in range(self.lmax+1):
            data = self._grep(tag="!      %d" % l).data
            assert l == int(data[0, 0])
            ae_atan_logder_l[l] = atan_logder(energies=data[:, 1], values=data[:, 2])
            ps_atan_logder_l[l] = atan_logder(energies=data[:, 1], values=data[:, 3])

        return self.AePsNamedTuple(ae=ae_atan_logder_l, ps=ps_atan_logder_l)

    @property
    def ene_vs_ecut(self):
        """Convergence of energy versus ecut for different l-values."""
        #convergence profiles, (ll=0,lmax)
        #!C     0    5.019345    0.010000
        #...
        #!C     1   19.469226    0.010000
        conv_data = collections.namedtuple("ConvData", "energies values")
        conv_l = {}

        for l in range(self.lmax+1):
            data = self._grep(tag="!C     %d" % l).data
            conv_l[l] = conv_data(energies=data[:, 1], values=data[:, 2])

        return conv_l

    def make_plotter(self, plotter_class=None):
        """
        Builds an instance of PseudoGenDataPlotter. One can customize the behavior
        of the plotter by passing a subclass via plotter_class.
        """
        plotter_class = PseudoGenDataPlotter if plotter_class is None else plotter_class
        kwargs = {k: getattr(self, k) for k in plotter_class.all_keys}

        return plotter_class(**kwargs)

    def _grep(self, tag, beg=0):
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
                # Exit because we know there's only one section starting with tag'
                if intag != -1:
                    stop = beg + i
                    break

        if not data:
            return self.GrepResults(data=None, start=intag, stop=stop)
        else:
            return self.GrepResults(data=np.array(data), start=intag, stop=stop)


def add_mpl_kwargs(method):
    """
    Decorate bound methods adding boilerplate code
    needed to customize matplotlib plots
    """
    from functools import wraps

    @wraps(method)
    def wrapper(*args, **kwargs):
        self, ax = args[0:2]

        #title = kwargs.pop("title", None)
        #if title is not None:
        #    ax.set_title(title)

        # Set yticks and labels.
        #ylabel = kwargs.pop("ylabel", None)
        #if ylabel is not None:
        #    ax.set_ylabel(ylabel)

        #xlabel = kwargs.pop("xlabel", None)
        #if xlabel is not None:
        #    ax.set_xlabel("xlabel")

        ax.grid(True)
        return method(self, ax, **kwargs)

    return wrapper


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
    linewidth, markersize = 2, 1

    linestyle_aeps = dict(ae="solid", ps="dashed")
    color_l = {-1: "blue", 0: "red", 1: "green", 2: "yellow", 3: "magenta"}

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

    def plot_key(self, key, **kwargs):
        """Plot a singol quantity specified by key."""
        fig = self._mplt.figure()
        ax = fig.add_subplot(1, 1, 1)
        # key --> self.plot_key()
        getattr(self, "plot_" + key)(ax, **kwargs)
        self._mplt.show()

    def _wf_pltopts(self, l, aeps):
        """Plot options for wavefunctions."""
        return dict(
            color=self.color_l[l], linestyle=self.linestyle_aeps[aeps],
            linewidth=self.linewidth, markersize=self.markersize)

    @add_mpl_kwargs
    def plot_atan_logders(self, ax, **kwargs):
        """Plot arctan of logder on axis ax."""
        ae, ps = self.atan_logders.ae, self.atan_logders.ps

        lines, legends = [], []
        for l, ae_alog in ae.items():
            ps_alog = ps[l]

            ae_line, = ax.plot(ae_alog.energies, ae_alog.values, **self._wf_pltopts(l, "ae"))
            ps_line, = ax.plot(ps_alog.energies, ps_alog.values, **self._wf_pltopts(l, "ps"))

            lines.extend([ae_line, ps_line])
            legends.extend(["AE l=%s" % str(l), "PS l=%s" % str(l)])

        ax.set_xlabel("Energy [Ha]")
        ax.set_title("ARCTAN(Log Derivatives)")
        ax.legend(lines, legends, loc="best", shadow=True)
        return lines

    @add_mpl_kwargs
    def plot_radial_wfs(self, ax, **kwargs):
        """
        Plot ae and ps radial wavefunctions on axis ax.

        lselect:
            List to select l channels
        """
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

        ax.set_xlabel("r [Bohr]")
        ax.set_title("Wave Functions")
        ax.legend(lines, legends, loc="best", shadow=True)
        return lines

    @add_mpl_kwargs
    def plot_projectors(self, ax, **kwargs):
        """
        Plot oncvpsp projectors on axis ax.

        lselect:
            List to select l channels
        """
        lselect = kwargs.get("lselect", [])

        lines, legends = [], []
        for nl, proj in self.projectors.items():
            if nl.l in lselect: continue
            line, = ax.plot(proj.rmesh, proj.values,
                            linewidth=self.linewidth, markersize=self.markersize)

            lines.append(line)
            legends.append("Proj %s" % str(nl))

        ax.set_xlabel("r [Bohr]")
        ax.set_title("Projector Wave Functions")
        ax.legend(lines, legends, loc="best", shadow=True)
        return lines

    @add_mpl_kwargs
    def plot_densities(self, ax, **kwargs):
        """Plot ae, ps and model densities on axis ax."""
        lines, legends = [], []
        for name, rho in self.densities.items():
            line, = ax.plot(rho.rmesh, rho.values,
                            linewidth=self.linewidth, markersize=self.markersize)

            lines.append(line)
            legends.append(name)

        ax.legend(lines, legends, loc="best", shadow=True)
        ax.set_xlabel("r [Bohr]")
        ax.set_title("Charge densities")
        return lines

    @add_mpl_kwargs
    def plot_potentials(self, ax, **kwargs):
        """Plot vl and vloc potentials on axis ax"""
        lines, legends = [], []
        for l, pot in self.potentials.items():
            line, = ax.plot(pot.rmesh, pot.values, **self._wf_pltopts(l, "ae"))
            lines.append(line)

            if l == -1:
                legends.append("Vloc")
            else:
                legends.append("PS l=%s" % str(l))

        ax.set_xlabel("r [Bohr]")
        ax.set_title("Ion Pseudopotentials")
        ax.legend(lines, legends, loc="best", shadow=True)
        return lines

    @add_mpl_kwargs
    def plot_ene_vs_ecut(self, ax, **kwargs):
        """Plot the converge of ene wrt ecut on axis ax."""
        lines, legends = [], []
        for l, data in self.ene_vs_ecut.items():
            line, = ax.plot(data.energies, data.values, **self._wf_pltopts(l, "ae"))

            lines.append(line)
            legends.append("Conv l=%s" % str(l))

        ax.set_xlabel("Ecut [Ha]")
        ax.set_title("Energy Error per Electron (Ha)")

        ax.legend(lines, legends, loc="best", shadow=True)
        ax.set_yscale("log")
        return lines

    def _finalize_fig(self, fig, **kwargs):
        """Boilerplate code performed before returning the matplotlib figure."""
        title = kwargs.get("title", None)
        if title is not None:
            fig.suptitle(title)

        if kwargs.get("show", True):
            self._mplt.show()

        savefig = kwargs.get("savefig", None)
        if savefig is not None:
            fig.savefig(savefig)

        return fig

    def plot_atanlogder_econv(self, **kwargs):
        """Plot atan(logder) and ecut converge on the same figure. Returns matplotlib Figure"""
        fig, ax_list = self._mplt.subplots(nrows=2, ncols=1, sharex=False, squeeze=False)
        ax_list = ax_list.ravel()

        self.plot_atan_logders(ax_list[0])
        self.plot_ene_vs_ecut(ax_list[1])

        return self._finalize_fig(fig, **kwargs)

    def plot_dens_and_pots(self, **kwargs):
        """Plot densities and potentials on the same figure. Returns matplotlib Figure"""
        fig, ax_list = self._mplt.subplots(nrows=2, ncols=1, sharex=False, squeeze=False)
        ax_list = ax_list.ravel()

        self.plot_densities(ax_list[0])
        self.plot_potentials(ax_list[1])

        return self._finalize_fig(fig, **kwargs)

    def plot_waves_and_projs(self, **kwargs):
        """Plot ae-ps wavefunctions and projectors on the same figure. Returns matplotlib Figure"""
        lmax = max(nl.l for nl in self.radial_wfs.ae.keys())
        fig, ax_list = self._mplt.subplots(nrows=lmax+1, ncols=2, sharex=True, squeeze=False)

        for l in range(lmax+1):
            ax_idx = lmax - l
            self.plot_radial_wfs(ax_list[ax_idx][0], lselect=[l])
            self.plot_projectors(ax_list[ax_idx][1], lselect=[l])

        return self._finalize_fig(fig, **kwargs)
