# coding: utf-8
"""Public API to access the pseudopotential tables."""
import os
import json
import tempfile
import logging
import numpy as np

from collections import OrderedDict
from monty.collections import AttrDict
from monty.string import list_strings
from monty.fnmatch import WildCard
from monty.termcolor import cprint
from monty.os.path import which
from pymatgen.core.periodic_table import Element
from pymatgen.core.xcfunc import XcFunc
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt
from abipy.flowtk.pseudos import Pseudo, PseudoTable
from pseudo_dojo.core.dojoreport import DojoReport


logger = logging.getLogger(__name__)


def dojopseudo_from_file(filepath):
    """
    Factory function used to construct a :class:`Pseudo` object from file.
    A DojoPseudo has a DojoReport section and this function adds the report
    to the object.

    Args:
        filepath: Path of the pseudopotential file or djrepo file.

    .. note::

        We cannot subclass Pseudo because it's actually the abstract base class
        and Pseudo.from_file is the factory function that returns the concreate subclass.
    """
    filepath = os.path.abspath(filepath)

    dojo_report = None
    if filepath.endswith(".djrepo"):
        dojo_report = DojoReport.from_file(filepath)
        pp_basename = dojo_report["basename"]
        filepath = os.path.join(os.path.dirname(filepath), pp_basename)

    # Init pseudo from file. Return None if parser error.
    pseudo = Pseudo.from_file(filepath)
    if pseudo is None: return pseudo
    #pseudo.__class__.dojo_report = property(lambda self: self.a + 1)

    if dojo_report is not None:
        # We've already read the report.
        pseudo.dojo_report = dojo_report
        return pseudo

    # Read DojoReport and add it to pseudos
    root, ext = os.path.splitext(filepath)
    djrepo = root + ".djrepo"
    if not os.path.exists(djrepo):
        raise RuntimeError("Cannot find djrepo file at %s" % djrepo)
    pseudo.dojo_report = DojoReport.from_file(djrepo)

    return pseudo


def add_dojoreport_to_pseudo(pseudo):
    """Add the dojoreport to a pseudo."""
    if not hasattr(pseudo, "dojo_report"):
        pseudo.dojo_report = DojoReport.from_file(pseudo.djrepo_path)
    return pseudo


class DojoInfo(AttrDict):
    """
    Dictionary with metadata associated to the PseudoDojo table.
    """
    # See https://github.com/Julian/jsonschema
    JSON_SCHEMA = {
        "$schema": "http://json-schema.org/schema#",

        "type": "object",
        "properties": {
            "pp_type": {"type": "string", "enum": ["NC", "PAW"]},
            "xc_name": {"type": "string", "enum": XcFunc.aliases()},
            "authors": {"type": "array"},
            "description": {"type": "string"},
            "references": {"type": "array"},
            "dojo_dir": {"type": "string"},
            #"generation_date": {"type": "string", "format": "date"},
            #"tags": {"type": "array", "items": {"type": "string", "enum": ["accuracy", "efficiency"]}},
            #"relativity": {"type": "string", "enum": ["non-relativistic", "scalar-relativistic", "relativistic"]}
        },
        "required": ["pp_type", "xc_name", "authors", "description", "references", "dojo_dir",
                    #"generation_date", "tags", "relativity"
                    ],
    }

    def __init__(self, *args, **kwargs):
        super(DojoInfo, self).__init__(*args, **kwargs)
        self["xc"] = XcFunc.from_name(self["xc_name"])

    def validate_json_schema(self):
        """Validate DojoInfo with validictory."""
        from jsonschema import validate
        validate(self, self.JSON_SCHEMA)

    @classmethod
    def get_template_dict(cls):
        """Return a dictionary with the keys that must be filled by the user."""
        return {k: str(v) for k, v in cls.JSON_SCHEMA["properties"].items()}

    @property
    def isnc(self):
        """True if norm-conserving pseudopotential."""
        return self.pp_type == "NC"

    @property
    def ispaw(self):
        """True if PAW pseudopotential."""
        return self.pp_type == "PAW"


class DojoTable(PseudoTable):
    """
    This a base class that is not supposed to be exposed in the public API.
    End-users will mainly use `OfficialDojoTable` istances.

    A DojoTable is a subclass of `PseudoTable` so that we can easily add extra properties or methods.
    needed by the pseudodojo. No restriction is enforced on the content of the table e.g. DojoTable
    can contain multiple pseudos for a given element, pseudos with mixed XC etc.
    All the methods implemented here must take this possibilty into account.
    Methods and properties that assume a well defined set of pseudos fulfilling the pseudo_dojo
    constraints should be implemented in OfficialDojoTable (see below).
    """

    def __init__(self, pseudos):
        super(DojoTable, self).__init__(pseudos)
        # Add dojo_report to pseudos.
        for p in self:
            add_dojoreport_to_pseudo(p)

    @classmethod
    def from_dojodir(cls, top, exclude_wildcard=None, exclude_basenames=None):
        """
        Initialize the table from one of the top level directories located
        in the pseudo_dojo.pseudos directory.

        Args:
            top: top level directory
            exclude_basenames: Optional string or list of strings with the
                pseudo basenames to be excluded.
            exclude_wildcard: String of tokens separated by "|". Each token represents a pattern.
                to be exluded
                Example:
                  wildcard="*_r.psp8|*.xml" selects only those files that do not end with _r.psp8 or .xml

        .. warning::

            The table may contain multiple pseudos for a given chemical element.
            Don't use this method unless you need this feature and you know what you are doing.
        """
        # Read metadata from the __init__.py file
        import imp
        module_name = os.path.join(top, "__init__.py")
        if not os.path.isfile(module_name):
            raise RuntimeError("__init_.py file is missing in dir: %s" % top)

        meta = imp.load_source(module_name, os.path.join(top, "__init__.py"))

        # Gather all pseudos starting from the current working directory
        all_symbols = set(e.symbol for e in Element)
        dirs = [os.path.join(top, d) for d in os.listdir(top) if d in all_symbols]

        exclude = set(list_strings(exclude_basenames)) if exclude_basenames is not None else set()

        paths = []
        for dr in dirs:
            paths.extend(os.path.join(dr, f) for f in os.listdir(dr)
                         if f.endswith(meta.pseudo_ext) and f not in exclude)

        if exclude_wildcard is not None:
            wild = WildCard(exclude_wildcard)
            paths = [p for p in paths if not wild.match(os.path.basename(p))]

        pseudos = []
        for p in paths:
            pseudo = dojopseudo_from_file(p)
            if pseudo is None:
                print("Error while parsing:", p)
                continue
            pseudos.append(pseudo)

        return cls(pseudos).sort_by_z()

    @classmethod
    def from_txtfile(cls, path):
        """
        Initialize the table from a text file containing the relative path
        of the pseudopotential (one name per line).
        """
        path = os.path.abspath(path)
        root, _ = os.path.split(path)
        pseudos = []
        with open(path, "rt") as fh:
            for line in fh:
                line = line.strip()
                if line.startswith("#"): continue
                p = os.path.join(root, line)
                pseudo = dojopseudo_from_file(p)
                if pseudo is None:
                    raise RuntimeError("Error while parsing:", p)
                pseudos.append(pseudo)

        return cls(pseudos).sort_by_z()

    def to_djson(self, verbose=0, ignore_dup=False):
        """
        Tool used by the PseudoDojo maintainers to build a dictionary
        with **partial** information on the table. This dictionary can be used as
        an initial template for the creation of a new djson file.

        Args:
            verbose: Verbosity level.
            ignore_dup: if set to True, duplicated pseudos are ignored, only
                    the first pseudo is reported. Use it wisely!
        """
        # Add template for dojo_info section
        d = {"dojo_info": DojoInfo.get_template_dict()}

        def djson_entry(p):
            jdict = p.as_dict()
            entry = OrderedDict([(k, jdict[k]) for k in ("basename", "Z_val", "l_max", "md5")])

            res = p.dojo_report.get_last_df_results(with_soc=False)
            entry["dfact_meV"] = res.get("dfact_meV", None)
            entry["dfactprime_meV"] = res.get("dfactprime_meV", None)
            entry["tags"] = p.dojo_report.get("tags", [])
            entry["hints"] = p.dojo_report.get("hints", {})

            return entry

        # Add pseudo_metadata section (sorted by Z).
        # If there are multiple pseudos per element, we create a list of dicts.
        # Authors of the table, will have to select one.
        d["pseudos_metadata"] = meta = OrderedDict()
        for p in self.sort_by_z():
            if verbose: print("Analyzing ",p)
            if p.symbol in meta:
                # Handle multiple pseudos.
                # THIS IS A HACK FOR GENERATING djson files for testing purpose.
                if ignore_dup: continue
                old = meta[p.symbol]
                if not isinstance(old, list): old = [old]
                old.append(djson_entry(p))
                meta[p.symbol] = old
            else:
                meta[p.symbol] = djson_entry(p)

        return d

    # TODO: Move to require_hints == True
    def dojo_find_errors(self, md5dict, require_hints=False):
        """
        This function tests whether the table fulfill the requirements
        imposed by the PseudoDojo. More specifically:

            #. One pseudo per element.

            #. All pseudos should have a valid dojo report with hints

            #. The md5 value computed from the pseudo potential file must agree
               with the one found in the djson file.

        Args:
            md5dict: dict basename --> md5 hash value

        Return:
            List of errors
        """
        errors = []
        eapp = errors.append

        # One pseudo per element.
        unique_symbols = set(p.symbol for p in self)
        if len(unique_symbols) != len(self):
            eapp("Found multiple pseudos for a given element.")

        # Test pseudopotential and dojo_report.
        for p in self:
            if not p.has_dojo_report:
                eapp("[%s] does not have the DOJO_REPORT section" % repr(p))
                continue

            estring = p.dojo_report.check()
            if estring:
                eapp("[%s]:" % repr(p) + estring)

            if require_hints and not p.dojo_report.has_hints:
                eapp("[%s] does not have hints" % repr(p))

            if p.md5 != md5dict[p.basename]:
                eapp("[%s] p.mdf5 != mdf5dict[p.basename]\n%s, %s" % (repr(p), p.md5, md5dict[p.basename]))

        # Test support for SOC. All the pseudos much have the same level.
        # At present, this check makes sense only for NC pseudos.
        # PAW pseudos support SOC within the on-site approach.
        for i, p in enumerate(self):
            if i == 0: p0 = p
            if p.supports_soc == p0.supports_soc: continue
            eapp("[%s] has different SOC characteristics" % repr(p))

        return errors

    def get_dojo_dataframe(self, **kwargs):
        """
        Buid a pandas :class:`DojoDataFrame` with the most important parameters extracted from the
        `DOJO_REPORT` section of each pseudo in the table.

        Returns:
            (frame, errors)

            where frame is the pandas :class:`DataFrame` and errors is a list of errors
            encountered while trying to read the `DOJO_REPORT` from the pseudopotential file.
        """
        from pseudo_dojo.core.dojoreport import DojoDataFrame
        return DojoDataFrame.from_pseudos(self)

    def get_dfgbrv_dataframe(self, raise_if_none_dojoreport=False):
        """
        Build and return a pandas :class:`DataFrame` in the form.

            basename     deltafactor  df_prime  gbrv_bcc  gbrv_fcc  symbol   md5
            H-high.psp8  0.074830     1.258340  0.028904  0.024726  H        5863396c90149cbe12af496141bde0d0
            ...

        where `gbrv_bcc` and `gbrv_fcc` are the relative errors (in percentage) wrt the AE calculations.

        Args:
            raise_if_none_dojoreport: If True, a ValueError is raised if one of the pseudo does not
                have the dojo_report else a warning is emitted.
        """
        from pseudo_dojo.core.dojoreport import DfGbrvDataFrame
        return DfGbrvDataFrame.from_pseudos(self, raise_if_none_dojoreport=raise_if_none_dojoreport)

    #def get_raren_dataframe(self):
    #    from pseudo_dojo.refdata.lantanides.database import raren_database
    #    xc = self
    #    for i, pseudo in enumerate(self):
    #        if i == 0:
    #            table = raren_database(pseudo.xc).table.copy()
    #        else
    #            assert table.xc == pseudo.xc

    #    return table

    def dojo_compare(self, what="all", **kwargs):
        """
        Compare ecut convergence and Deltafactor, GBRV results
        """
        import matplotlib.pyplot as plt
        show = kwargs.pop("show", True)
        what = list_strings(what)
        figs = []

        if all(p.dojo_report.has_trial("deltafactor") for p in self) and \
           any(k in what for k in ("all", "ecut")):

            fig_etotal, ax_list = plt.subplots(nrows=len(self), ncols=1, sharex=True, squeeze=True)
            figs.append(fig_etotal)

            for ax, pseudo in zip(ax_list, self):
                pseudo.dojo_report.plot_etotal_vs_ecut(ax=ax, show=False, label=pseudo.basename)
            if show: plt.show()

        if all(p.dojo_report.has_trial("deltafactor") for p in self) and \
           any(k in what for k in ("all", "df", "deltafactor")):

            fig_deltafactor, ax_grid = plt.subplots(nrows=5, ncols=len(self), sharex=True, sharey="row", squeeze=False)
            figs.append(fig_deltafactor)

            for ax_list, pseudo in zip(ax_grid.T, self):
                pseudo.dojo_report.plot_deltafactor_convergence(xc=pseudo.xc, ax_list=ax_list, show=False)

            fig_deltafactor.suptitle(" vs ".join(p.basename for p in self))
            if show: plt.show()

        # Compare GBRV results
        if all(p.dojo_report.has_trial("gbrv_bcc") for p in self) and \
           any(k in what for k in ("all", "gbrv")):

            fig_gbrv, ax_grid = plt.subplots(nrows=2, ncols=len(self), sharex=True, sharey="row", squeeze=False)
            figs.append(fig_gbrv)

            for ax_list, pseudo in zip(ax_grid.T, self):
                pseudo.dojo_report.plot_gbrv_convergence(ax_list=ax_list, show=False)

            fig_gbrv.suptitle(" vs ".join(p.basename for p in self))
            if show: plt.show()

        return figs

    def dojo_nbcompare(self, what="all", **kwargs):
        """
        Generate an ipython notebook to compare the results in the dojoreport (calls dojo_compare).
        """
        paths = [p.path for p in self]
        import nbformat
        nbf = nbformat.v4
        nb = nbf.new_notebook()

        nb.cells.extend([
            #nbf.new_markdown_cell("# This is an auto-generated notebook for %s" % os.path.basename(pseudopath)),
            nbf.new_code_cell("""\
from __future__ import print_function, division, unicode_literals
%matplotlib notebook"""),

            nbf.new_code_cell("""\
from pseudo_dojo.core.pseudos import DojoTable
pseudos = DojoTable(%s)""" % str(paths)),
            nbf.new_code_cell("pseudos.dojo_compare(what='%s')" % what),
        ])

        _, nbpath = tempfile.mkstemp(suffix='.ipynb', text=True)

        import io
        with io.open(nbpath, 'wt', encoding="utf8") as f:
            nbformat.write(nb, f)

        if which("jupyter") is None:
            raise RuntimeError("Cannot find jupyter in PATH. Install it with `pip install`")
        return os.system("jupyter notebook %s" % nbpath)

    @add_fig_kwargs
    def plot_dfgbrv_dist(self, **kwargs):
        """
        Plot four distribution plots for the deltafactor, deltafactor prime and the
        relative errors for the GBRV fcc, bcc structures.

        Return: `matplotlib` figure.
        """
        frame = self.get_dfgbrv_dataframe()
        return frame.plot_dfgbrv_dist(**kwargs)

    @add_fig_kwargs
    def plot_scalar_vs_fully_relativistic(self, what="df", **kwargs):
        # Build pandas dataframe with results.
        rows = []
        for p in self:
            if not p.has_dojo_report:
                cprint("Cannot find dojo_report in %s" % p.basename, "magenta")
                continue
            report = p.dojo_report
            row = {att: getattr(p, att) for att in ("basename", "Z", "Z_val", "l_max")}

            # Get deltafactor data with/without SOC
            soc_dict = p.dojo_report.get_last_df_results(with_soc=True)
            sr_dict = p.dojo_report.get_last_df_results(with_soc=False)
            row.update(sr_dict)
            row.update(soc_dict)

            # Get GBRV results for bcc and fcc.
            for struct_type in ("bcc", "fcc"):
                soc_dict = p.dojo_report.get_last_gbrv_results(struct_type, with_soc=True)
                sr_dict = p.dojo_report.get_last_gbrv_results(struct_type, with_soc=False)
                row.update(sr_dict)
                row.update(soc_dict)

            rows.append(row)

        import pandas as pd
        frame = pd.DataFrame(rows)

        def print_frame(x):
            with pd.option_context('display.max_rows', len(x),
                                   'display.max_columns', len(list(x.keys()))):
                print(x)

        # Create axes
        import matplotlib.pyplot as plt

        if what == "df":
            fig, ax_list = plt.subplots(nrows=2, ncols=2, sharex=True, squeeze=False)
            ax_list = ax_list.ravel()

            # Compute absolute differences SOC - SR
            for i, vname in enumerate(["dfact_meV", "dfactprime_meV",
                                      "v0", "b0_GPa"]):
                newcol = "dsoc_" + vname
                frame[newcol] = frame[vname + "_soc"] - frame[vname]
                print_frame(frame[["basename", "Z", newcol]])
                frame.plot(x="Z", y=newcol, ax=ax_list[i], kind="scatter", grid=True)
                #frame.plot.scatter(x="Z", y=newcol, s=20*(frame["dfact_meV"] + 1),
                #                    ax=ax_list[i], grid=True)

            #print("entries without V0")
            #print_frame(frame[frame["v0"].isnull()])
            #print("entries without V0_soc")
            #print_frame(frame[frame["v0_soc"].isnull()])

        elif what == "gbrv":
            # Plot GBRV results with/without SOC
            fig, ax_list = plt.subplots(nrows=2, ncols=1, sharex=True, squeeze=False)
            ax_list = ax_list.ravel()

            """
            for i, stype in enumerate(("bcc", "fcc")):
                # "gbrv_fcc_soc_a0_rel_err"
                for soc, color in zip(("", "_soc"), ("blue", "red")):
                    col = "gbrv_" + stype + soc + "_a0_rel_err"
                    #print_frame(frame)
                    frame.plot(x="Z", y=col, ax=ax_list[i], kind="scatter", grid=True, color=color)
            """

            # Plot GBRV diff between SR and FR with SOC
            #frame["dsoc_bcc"] = frame["gbrv_bcc_soc_a0_rel_err"] - frame["gbrv_bcc_a0_rel_err"]
            #frame.plot(x="Z", y="dsoc_bcc", ax=ax_list[0], kind="scatter", grid=True)
            #frame["dsoc_fcc"] = frame["gbrv_fcc_soc_a0_rel_err"] - frame["gbrv_fcc_a0_rel_err"]
            #frame.plot(x="Z", y="dsoc_fcc", ax=ax_list[1], kind="scatter", grid=True)

            frame["dsoc_bcc"] = \
                100 * (frame["gbrv_bcc_soc_a0"] - frame["gbrv_bcc_a0"]) / frame["gbrv_bcc_a0"]
            frame.plot(x="Z", y="dsoc_bcc", ax=ax_list[0], kind="scatter", grid=True)
            frame["dsoc_fcc"] = \
                100 * (frame["gbrv_fcc_soc_a0"] - frame["gbrv_fcc_a0"]) / frame["gbrv_fcc_a0"]
            frame.plot(x="Z", y="dsoc_fcc", ax=ax_list[1], kind="scatter", grid=True)

            #print_frame(frame[frame["gbrv_fcc_soc_a0_rel_err"].isnull()])
            #print_frame(frame[frame["gbrv_fcc_a0_rel_err"].isnull()])
            #print_frame(frame[["basename", "Z", "Z_val", "dsoc_bcc", "dsoc_fcc"]])
            #print_frame(frame)
        else:
            raise ValueError("Unsupported option %s" % what)

        for ax in ax_list: ax.set_xlim(0)

        return fig

    @add_fig_kwargs
    def plot_hints(self, with_soc=False, **kwargs):
        # Build pandas dataframe with results.
        rows = []
        for p in self:
            if not p.has_dojo_report:
                cprint("Cannot find dojo_report in %s" % p.basename, "magenta")
                continue
            report = p.dojo_report
            row = {att: getattr(p, att) for att in ("basename", "symbol", "Z", "Z_val", "l_max")}

            # Get deltafactor data with/without SOC
            df_dict = report.get_last_df_results(with_soc=with_soc)
            row.update(df_dict)
            for struct_type in ["fcc", "bcc"]:
                gbrv_dict = report.get_last_gbrv_results(struct_type, with_soc=with_soc)
            row.update(gbrv_dict)

            # Get the hints
            hint = p.hint_for_accuracy(accuracy="normal")
            row.update(dict(ecut=hint.ecut, pawecutdg=hint.pawecutdg))

            rows.append(row)

        import pandas as pd
        frame = pd.DataFrame(rows)

        def print_frame(x):
            import pandas as pd
            with pd.option_context('display.max_rows', len(x),
                                   'display.max_columns', len(list(x.keys()))):
                print(x)

        print_frame(frame)
        # Create axes
        #import matplotlib.pyplot as plt

        import seaborn as sns
        ax, fig, plt = get_ax_fig_plt(ax=None)

        #order = sort_symbols_by_Z(set(frame["element"]))

        # Box plot
        ax = sns.boxplot(x="symbol", y="ecut", data=frame, ax=ax, #order=order,
                         whis=np.inf, color="c")
        # Add in points to show each observation
        sns.stripplot(x="symbol", y="ecut", data=frame, ax=ax, #order=order,
                      jitter=True, size=5, color=".3", linewidth=0)

        sns.despine(left=True)
        ax.set_ylabel("Relative error %")
        ax.grid(True)

        return fig

    def make_open_notebook(self, nbpath=None, foreground=False):
        """
        Generate an ipython notebook and open it in the browser.

        Args:
            nbpath: If nbpath is None, a temporay file is created.
            foreground: By default, jupyter is executed in background and stdout, stderr are redirected
            to devnull. Use foreground to run the process in foreground

        Return:
            system exit code.

        Raise:
            RuntimeError if jupyter is not in $PATH
        """
        nbpath = self.write_notebook(nbpath=nbpath)

        if which("jupyter") is None:
            raise RuntimeError("Cannot find jupyter in PATH. Install it with `pip install`")

        if foreground:
            cmd = "jupyter notebook %s" % nbpath
            return os.system(cmd)

        else:
            cmd = "jupyter notebook %s &> /dev/null &" % nbpath
            print("Executing:", cmd)

            import subprocess
            try:
                from subprocess import DEVNULL # py3k
            except ImportError:
                DEVNULL = open(os.devnull, "wb")

            process = subprocess.Popen(cmd.split(), shell=False, stdout=DEVNULL) #, stderr=DEVNULL)
            cprint("pid: %s" % str(process.pid), "yellow")

    def write_notebook(self, nbpath=None):
        """
        Write an ipython notebook to nbpath
        If nbpath is None, a temporay file is created.
        Return path to the notebook.

        See also:
            http://nbviewer.jupyter.org/github/maxalbert/auto-exec-notebook/blob/master/how-to-programmatically-generate-and-execute-an-ipython-notebook.ipynb
        """
        # Get frame and write data in JSON format to tmp file so that we can reread in the notebook.
        frame, errors = self.get_dojo_dataframe()

        _, json_path = tempfile.mkstemp(suffix='.json', text=True)
        frame.to_json(path_or_buf=json_path)

        import nbformat
        nbf = nbformat.v4
        nb = nbf.new_notebook()

        nb.cells.extend([
            nbf.new_markdown_cell("# This is an auto-generated notebook"),
            nbf.new_code_cell("""\
from __future__ import print_function, division, unicode_literals

%matplotlib notebook
import seaborn as sns
from tabulate import tabulate
from IPython.display import display
from pseudo_dojo.core.dojoreport import DojoDataFrame"""),

            nbf.new_markdown_cell("## Init table from filenames:"),
            nbf.new_code_cell("dojo_frame = DojoDataFrame.from_json_file('%s')" % json_path),
            nbf.new_code_cell("#display(best_frame)"),

            nbf.new_code_cell("best_frame = dojo_frame.select_best()"),
            nbf.new_code_cell("display(best_frame)"),
            nbf.new_code_cell("tabulate(best_frame.describe(), headers='keys')"),

            nbf.new_code_cell("fig = options.pseudos.plot_dfgbrv_dist()"),

            nbf.new_code_cell("""\
for row in dojo_frame.myrows():
    print("row:", row)
    row_frame = dojo_frame.select_rows(row)
    row_frame.plot_hist()"""),

            nbf.new_code_cell("""\
for family in dojo_frame.myfamilies():
    print("family:", family)
    family_frame = dojo_frame.select_family(family)
    family_frame.plot_hist()"""),
        ])

        if nbpath is None:
            _, nbpath = tempfile.mkstemp(suffix='.ipynb', text=True)

        import io
        with io.open(nbpath, 'wt', encoding="utf8") as f:
            nbformat.write(nb, f)

        return nbpath


class OfficialDojoTable(DojoTable):
    """
    An official pseudopotential table provided by the pseudo_dojo.
    The table contains **one** pseudo per element and can be used for production.
    Client code usually access these tables via the OfficialTables() interface.
    Note the followig constraints enforced in OfficialDojoTable:

    #. One pseudo per element.
        #. Each pseudo has been generated with the same XC and with the same
           treatment of relativitistic corrections.
        #. All pseudos must have a valid dojo report with hints
        #. The md5 value computed from the pseudo potential file must agree
           with the one found in the djson file.

    .. attribute:: xc

        XC functional
    """
    @classmethod
    def from_dojodir(cls,dojodir,accuracy='standard'):
        """Use a dojodir string to get a djson file and initialize the class"""
        import glob
        from pseudo_dojo.pseudos import dojotable_absdir
        dojodir_path = dojotable_absdir(dojodir)
        djson_path = os.path.join(dojodir_path,accuracy+'.djson')
        if not os.path.isfile(djson_path):
            filenames = glob.glob(os.path.join(dojodir_path,"*.djson"))
            accuracies = [os.path.basename(filename).replace('.djson','') for filename in filenames]
            raise FileNotFoundError("File {} does not exist. "
                                    "Found djson files for accuracy = {}".format(djson_path,accuracies))
        return cls.from_djson_file(djson_path)

    @classmethod
    def from_djson_file(cls, json_path):
        """
        Initialize the pseudopotential table from one of **official** djson files
        located in one of the subdirectories inside pseudo_dojo.pseudos.

        json_path contains the following dictionary in JSON format:

        {
        "dojo_info": {
              "pp_type": "NC",
              "xc_name": "PBE",
              "authors": ["J. Doe",],
              "generation_date": "2015-07-20",
              "description": "String",
              "tags": ["accuracy", "tag2"],
              "reference": "paper",
              "dojo_dir": "ONCVPSP-PBE",
        },
        "pseudos_metadata": {
            "Si": {
                "basename": "Si-dloc.psp8",
                "Z_val": 4.0,
                "l_max": 2,
                "md5": "ececcf5b26f34676694b630d6bc809ff"
            },
            "O": {
                "basename": "O-dmax.psp8",
                "Z_val": 6.0,
                "l_max": 2,
                "md5": "f7d0f3573362d89c81c41fc6b7b3e6ab"
            }
        }
        }
        """
        json_path = os.path.abspath(json_path)
        with open(json_path, "rt") as fh:
            d = json.load(fh)

        # Read and validate dojo_info.
        dojo_info = DojoInfo(**d["dojo_info"])
        try:
            dojo_info.validate_json_schema()
        except Exception as exc:
            print("Validation error in %s" % json_path)
            raise exc

        meta = d["pseudos_metadata"]

        if dojo_info.get("dojo_dir", None):
            from pseudo_dojo.pseudos import dojotable_absdir
            top = dojotable_absdir(dojo_info.dojo_dir)
        else:
            top = os.path.dirname(json_path)
        paths, md5dict = [], {}
        for esymb, m in meta.items():
            if isinstance(m, (list, tuple)):
                raise TypeError("Invalid djson file. Expecting dict but got list (multiple pseudos):\n\n %s" % str(m))

            path = os.path.join(top, esymb, m["basename"])
            paths.append(path)
            md5dict[m["basename"]] = m["md5"]

        # TODO: Avoid parsing the pseudos. Construct them from dict.
        new = cls(paths).sort_by_z()
        new.set_dojo_info(dojo_info)

        # TODO: To be activated
        #errors = new.dojo_find_errors(md5dict)
        #if errors:
        #    raise ValueError("\n".join(errors))

        return new

    @property
    def xc(self):
        """
        The :class:`XcFunc` object describing the XC functional used to generate the table.
        """
        return self._dojo_info.xc

    @property
    def dojo_info(self):
        try:
            return self._dojo_info
        except AttributeError:
            return {}

    #@dojo_info.setter
    def set_dojo_info(self, dojo_info):
        self._dojo_info = dojo_info
