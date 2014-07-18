"""Classes and functions for the post-processing of the results produced by ONCVPSP"""
from __future__ import print_function, division

import os
import collections
import numpy as np

from pprint import pprint
from pseudo_dojo.core import RadialFunction, RadialWaveFunction


class OncvOutputParserError(Exception):
    """Exceptions raised by OncvParser."""


class OncvOuptputParser(object):
    """Object to read and extract data from the output file of oncvpsp."""
    Error = OncvOutputParserError

    GrepResults = collections.namedtuple("GrepResults", "data, start, stop")

    # Used to store ae-pp wavefunctions, ae-pp log-derivatives
    AePsNamedTuple = collections.namedtuple("AePsNamedTuple", "ae, ps")

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

        # Read lmax (not very robust because we assume the user didn't change the template but oh well)
        for i, line in enumerate(self.lines):
            if line.startswith("# lmax"):
                self.lmax = int(self.lines[i+1])
                break
        else:
            raise self.Error("Cannot find #lmax line in output file %s" % self.filepath)

    def __str__(self):
        lines = []
        app = lines.append
        app("%s\nversion: %s\ndate: %s" % (self.calc_type, self.version, self.gendate))

        return "\n".join(lines)

    @property
    def potentials(self):
        """Radial functions with the non-local and local potentials."""
        #radii, charge, pseudopotentials (ll=0, 1, lmax)
        #!p   0.0099448   4.7237412  -7.4449470 -14.6551019
        vl_data = self._grep_nparr("!p").data
        lmax = len(vl_data[0]) - 3
        assert lmax == self.lmax

        # From 0 up to lmax
        ionpots_l = {}
        for l in range(lmax+1):
            ionpots_l[l] = RadialFunction("Ion Pseudopotential, l=%d" % l, vl_data[:, 0], vl_data[:, 2+l])

        # Local part is stored with l == -1
        vloc = self._grep_nparr("!L").data
        ionpots_l[-1] = RadialFunction("Local part, l=%d" % -1, vloc[:, 0], vloc[:, 1])

        return ionpots_l

    @property
    def densities(self):
        """Dictionary with charge densities on the radial mesh."""
        # radii, charge, core charge, model core charge
        # !r   0.0100642   4.7238866  53.4149287   0.0000000
        rho_data = self._grep_nparr("!r").data

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
            g = self._grep_nparr("&", beg=beg)
            #print(g)
            if g.data is None:
                break
            beg = g.stop + 1

            header = self.lines[g.start-2]
            n, l = header.split(",")[0:2]
            n = int(n.split("=")[1])
            l = int(l.split("=")[1])
            state = (n, l)
            print("Got state: %s" % str(state))

            rmesh = g.data[:, 1]
            ae_wf = g.data[:, 2]
            ps_wf = g.data[:, 3]

            ae_waves[state] = RadialWaveFunction(state, state, rmesh, ae_wf)
            ps_waves[state] = RadialWaveFunction(state, state, rmesh, ps_wf)

        return self.AePsNamedTuple(ae=ae_waves, ps=ps_waves)

    @property
    def projector_waves(self):
        """Read the radial wavefunctions."""
        # TODO: Check this
        #n= 1 2  l= 0, projecctor pseudo wave functions, well or 2nd valence
        #
        #@     0    0.009945    0.015274   -0.009284
        projectors_nl = {}
        beg = 0
        while True:
            g = self._grep_nparr("@", beg=beg)
            if g.data is None:
                break
            beg = g.stop + 1
            #header = self.lines[g.start-2]
            #n, l = header.split(",")[0:2]
            #n = int(n.split("=")[1])
            #l = int(l.split("=")[1])

            l = int(g.data[0, 0])
            rmesh = g.data[:, 1]

            #pprint(g.data)
            #print(range(len(g.data[0]) - 3))
            for n in range(len(g.data[0]) - 2):
                state = (n, l)
                print("Got projector with (n, l) state: %s" % str(state))
                projectors_nl[state] = RadialWaveFunction(state, state, rmesh, g.data[:, n+2])

        return projectors_nl

    @property
    def atan_logders(self):
        #log derivativve data for plotting, l= 0
        #atan(r * ((d psi(r)/dr)/psi(r))), r=  1.60
        #l, energy, all-electron, pseudopotential
        #
        #!      0    2.000000    0.706765    0.703758
        atan_logder = collections.namedtuple("AtanLogDer", "energies values")
        ae_atan_logder_l, ps_atan_logder_l = {}, {}

        for l in range(self.lmax+1):
            data = self._grep_nparr(tag="!      %d" % l).data
            assert l == int(data[0, 0])
            ae_atan_logder_l[l] = atan_logder(energies=data[:, 1], values=data[:, 2])
            ps_atan_logder_l[l] = atan_logder(energies=data[:, 1], values=data[:, 3])

        return self.AePsNamedTuple(ae=ae_atan_logder_l, ps=ps_atan_logder_l)

    @property
    def ene_vs_ecut(self):
        #convergence profiles, (ll=0,lmax)
        #!C     0    5.019345    0.010000
        #...
        #!C     1   19.469226    0.010000
        conv_data = collections.namedtuple("ConvData", "energies values")
        conv_l = {}

        for l in range(self.lmax+1):
            data = self._grep_nparr(tag="!C     %d" % l).data
            conv_l[l] = conv_data(energies=data[:, 1], values=data[:, 2])

        return conv_l

    def make_plotter(self, plotter_class=None):
        plotter_class = Plotter if plotter_class is None else plotter_class
        kwargs = {k: getattr(self, k) for k in plotter_class.all_keys}

        return plotter_class(**kwargs)

    def _grep_nparr(self, tag, beg=0):
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
            #raise ValueError("Cannot find tag %s in output file %s" % (tag, self.filepath))
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

        title = kwargs.pop("title", None)
        if title is not None:
            ax.set_title(title)

        # Set yticks and labels.
        ylabel = kwargs.pop("ylabel", None)
        if ylabel is not None:
            ax.set_ylabel(ylabel)

        xlabel = kwargs.pop("xlabel", None)
        if xlabel is not None:
            ax.set_xlabel("xlabel")

        ax.grid(True)

        return method(self, ax, **kwargs)

    return wrapper


class Plotter(object):

    all_keys = [
        "radial_wfs",
        #"projector_waves",
        "densities",
        "potentials",
        "atan_logders",
        "ene_vs_ecut",
    ]

    linestyle_aeps = dict(ae="solid", ps="dashed")
    color_l = {-1: "blue", 0: "red", 1: "green", 2: "yellow", 3: "magenta"}

    # Plot parameters
    linewidth = 2
    markersize = 1

    def __init__(self, **kwargs):
        import matplotlib.pyplot as _mplt
        self._mplt = _mplt

        for k in self.all_keys:
            setattr(self, k, kwargs.pop(k, {}))
        if kwargs:
            raise ValueError("Unknown keys: %s" % list(kwargs.keys()))

    def keys(self):
        return (k for k in self.all_keys if getattr(self, k))

    def plot_single(self, key, **kwargs):
        fig = self._mplt.figure()
        ax = fig.add_subplot(1, 1, 1)
        self._plot_driver(ax, key, **kwargs)
        self._mplt.show()

    def plot_slideshow(self, **kwargs):
        for key in self.keys():
            fig = self._mplt.figure()
            ax = fig.add_subplot(1, 1, 1)
            self._plot_driver(ax, key, **kwargs)
            self._mplt.show()

    def plot_all(self, **kwargs):
        num_plots = len(list(self.keys()))

        fig = self._mplt.figure()
        for ip, key in enumerate(self.keys()):
            ax = fig.add_subplot(num_plots, 1, ip)
            self._plot_driver(ax, key, **kwargs)

        self._mplt.show()

    def _plot_driver(self, ax, key, **kwargs):
        # key --> self.plot_key()
        return getattr(self, "plot_" + key)(ax, **kwargs)

    def _wf_pltopts(self, l, aeps):
        return dict(
            color=self.color_l[l], linestyle=self.linestyle_aeps[aeps],
            linewidth=self.linewidth, markersize=self.markersize)

    @add_mpl_kwargs
    def plot_atan_logders(self, ax, **kwargs):
        ae, ps = self.atan_logders.ae, self.atan_logders.ps

        lines, legends = [], []
        for l, ae_alog in ae.items():
            ps_alog = ps[l]

            ae_line, = ax.plot(ae_alog.energies, ae_alog.values, **self._wf_pltopts(l, "ae"))
            ps_line, = ax.plot(ps_alog.energies, ps_alog.values, **self._wf_pltopts(l, "ps"))

            lines.extend([ae_line, ps_line])
            legends.extend(["AE l=%s" % str(l), "PS l=%s" % str(l)])

        ax.legend(lines, legends, loc="best", shadow=True)
        return lines

    @add_mpl_kwargs
    def plot_radial_wfs(self, ax,  **kwargs):
        ae_wfs, ps_wfs = self.radial_wfs.ae, self.radial_wfs.ps

        lines, legends = [], []
        for state, ae_wf in ae_wfs.items():
            ps_wf = ps_wfs[state]
            # TODO Use namedtuple
            n, l = state[0], state[1]

            ae_line, = ax.plot(ae_wf.rmesh, ae_wf.values, **self._wf_pltopts(l, "ae"))
            ps_line, = ax.plot(ps_wf.rmesh, ps_wf.values, **self._wf_pltopts(l, "ps"))

            lines.extend([ae_line, ps_line])
            legends.extend(["AE l=%s" % str(l), "PS l=%s" % str(l)])

        ax.legend(lines, legends, loc="best", shadow=True)
        return lines

    @add_mpl_kwargs
    def plot_densities(self, ax,  **kwargs):

        lines, legends = [], []
        for name, rho in self.densities.items():
            line, = ax.plot(rho.rmesh, rho.values,
                            linewidth=self.linewidth, markersize=self.markersize)

            lines.append(line)
            legends.append(name)

        ax.legend(lines, legends, loc="best", shadow=True)
        return lines

    @add_mpl_kwargs
    def plot_potentials(self, ax,  **kwargs):

        lines, legends = [], []
        for l, pot in self.potentials.items():
            line, = ax.plot(pot.rmesh, pot.values, **self._wf_pltopts(l, "ae"))

            lines.append(line)
            legends.append("PS l=%s" % str(l))

        ax.legend(lines, legends, loc="best", shadow=True)
        return lines

    @add_mpl_kwargs
    def plot_ene_vs_ecut(self, ax, **kwargs):

        lines, legends = [], []
        for l, data in self.ene_vs_ecut.items():
            line, = ax.plot(data.energies, data.values, **self._wf_pltopts(l, "ae"))

            lines.append(line)
            legends.append("Conv l=%s" % str(l))

        ax.legend(lines, legends, loc="best", shadow=True)
        ax.set_yscale("log")
        return lines


def main():
    import argparse
    parser = argparse.ArgumentParser(epilog="",
                                     formatter_class=argparse.RawDescriptionHelpFormatter,)

    parser.add_argument('filename', default="", help="Path to the output file")

    parser.add_argument("-p", "--plot-mode", default="all",
                        help="Quantity to plot (default to all) Can be: %s" % str(["all", "slide"] + Plotter.all_keys))

    options = parser.parse_args()

    parser = OncvOuptputParser(options.filename)

    p = parser.make_plotter()

    if options.plot_mode == "all":
        p.plot_all()
    elif options.plot_mode == "slide":
        p.plot_slideshow()
    else:
        p.plot_single(key=options.plot_mode)


if __name__ == "__main__":
    main()
