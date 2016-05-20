# coding: utf-8
"""Public API to access the pseudopotential tables."""
from __future__ import division, print_function, unicode_literals

import os
import json
import numpy as np

from monty.collections import AttrDict
from monty.functools import lazy_property
from monty.string import list_strings
from pymatgen.core.periodic_table import PeriodicTable
#from pymatgen.core.xcfunc import XcFunc
from pymatgen.util.plotting_utils import add_fig_kwargs #, get_ax_fig_plt
from pymatgen.io.abinit.pseudos import PseudoTable


class DojoInfo(AttrDict):
    """
    Dictionary with metadata associated to the PseudoDojo table.
    """
    # See http://validictory.readthedocs.org/en/latest/usage.html#schema-options
    JSON_SCHEMA = {
        "type": "object",
        "properties": {
            "pseudo_type": {"type": "string", "enum": ["NC", "PAW"]},
            #"xc_name": {"type": "string", "enum": XcFunc.aliases()},
            "authors": {"type": "array"},
            #"generation_date": {"type": "string", "format": "date"},
            "description": {"type": "string"},
            "reference": {"type": "string"},
            "dojo_dir": {"type": "string"},
            #"tags": {"type": "array", "items": {"type": "string", "enum": ["accuracy", "efficiency"]}},
            #non-relativistic, scalar-relativistic or relativistic
        },
    }

    def validate_json_schema(self):
        """Validate DojoInfo with validictory."""
        import validictory
        validictory.validate(self, self.JSON_SCHEMA)

    @classmethod
    def get_template_dict(cls):
        """Return a dictionary with the keys that must be filled by the user."""
        return {k: str(v) for k, v in cls.JSON_SCHEMA["properties"].items()}

    @property
    def isnc(self):
        """True if norm-conserving pseudopotential."""
        return self.pseudo_type == "NC"

    @property
    def ispaw(self):
        """True if PAW pseudopotential."""
        return self.pseudo_type == "PAW"


class DojoTable(PseudoTable):
    """
    A pseudopotential table provided by the pseudo_dojo.
    We subclass `PseudoTable` so that we can easily add 
    extra properties or methods, if needed.
    """

    @classmethod
    def from_dojodir(cls, top, exclude_basenames=None):
        """
        Initialize the table of pseudos from one of the top level directories located
        in the pseudo_dojo.pseudos directory.

        Args:
            top: top level directory
            exclude_basenames: Optional string or list of strings with the
                pseudo basenames to be excluded.

        .. warning::

            The table may contain multiple pseudos for a given chemical element.
            Don't use this method unless you need this feature and you know what you are doing.
        """
        # Read metadata from the __init__.py file
        import imp
        module_name = os.path.join(top, "__init__.py")
        if not os.path.isfile(module_name):
            raise RuntimeError("__init_.py file is missing in dir: %s" % top)

        meta = imp.load_source(module_name, os.path.join(top, "__init__.py") )

        # Gather all pseudos starting from the current working directory
        all_symbols = set(element.symbol for element in PeriodicTable().all_elements)
        dirs = [os.path.join(top, d) for d in os.listdir(top) if d in all_symbols]

        exclude = set(list_strings(exclude_basenames)) if exclude_basenames is not None else set()

        paths = []
        for dr in dirs:
            paths.extend(os.path.join(dr, f) for f in os.listdir(dr)
                         if f.endswith(meta.pseudo_ext) and f not in exclude)

        return cls(paths).sort_by_z()

    @classmethod
    def from_djson_file(cls, djson_path):
        """
        Initialize the pseudopotential table from one of **official** djson files
        located in one of the subdirectories inside pseudo_dojo.pseudos.

        .. important::

            The table contains one pseudo per element and can be used for production calculations

        djson_path contains the dictionary:

        {
        "dojo_info": {
              "pseudo_type": "NC",
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
        with open(djson_path, "rt") as fh:
            d = json.load(fh)

        dojo_info = DojoInfo(**d["dojo_info"])
        dojo_info.validate_json_schema()
        #print(list(d.keys()))
        meta = d["pseudos_metadata"]

        top = os.path.dirname(djson_path)
        paths, md5dict = [], {}
        for esymb, m in meta.items():
            paths.append(os.path.join(top, esymb, m["basename"]))
            md5dict[m["basename"]] = m["md5"]

        new = cls(paths).sort_by_z()
        new.set_dojo_info(dojo_info)

        # TODO: To be activated
        #errors = new.dojo_find_errors(md5dict)
        #if errors:
        #    raise ValueError("\n".join(errors))

        return new

    def to_djson(self, **kwargs):
        """
        Build and return a dictionary with **partial** information
        on the table. This dictionary can be used as template for
        the creation of a new djson file.
        """
        # Add template for dojo_info section
        d = {"dojo_info": DojoInfo.get_template_dict()}

        def djson_entry(p):
            jdict = p.as_dict()
            return {k: jdict[k] for k in ["basename", "Z_val", "l_max", "md5"]}

        # Add pseudo_metadata section.
        # If there are multiple pseudos per element, we create a list of dicts.
        # Authors of the table, will select one.
        d["pseudos_metadata"] = meta = {}
        for p in self:
            if p.symbol in meta:
                continue # FIXME
                old = meta[p.symbol]
                if not isinstance(old, list): old = [old]
                old.append(djson_entry(p))
                meta[p.symbol] = old
            else:
                meta[p.symbol] = djson_entry(p)

        return d

    #@lazy_property
    #def xc(self):
    #    """The XcFunc object describing the XC functional used to generate the table."""
    #    xc = self[0].xc
    #    if any(p.xc != xc for p in self):
    #        raise TypeError("Found pseudos generated with different xc functionals")
    #    return xc

    @property
    def dojo_info(self):
        try:
            return self._dojo_info
        except AttributeError:
            return {}

    #@dojo_info.setter
    def set_dojo_info(self, dojo_info):
        self._dojo_info = dojo_info

    #def show_dojo_info(self):
    #    pprint(self.dojo_info)

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
                eapp("%s does not have the DOJO_REPORT section" % repr(p))
                continue

            estring = p.dojo_report.check()
            if estring:
                eapp(estring)

            if require_hints and not p.dojo_report.has_hints:
                eapp("%s does not have hints" % repr(p))

            if p.md5 != md5dict[p.basename]:
                eapp("p.mdf5 != mdf5dict[p.basename]\n%s, %s" % (p.md5, md5dict[p.basename]))

        return errors

    def get_dfgbrv_dataframe(self):
        """
        Build and return a pandas :class:`DataFrame` in the form.

            basename     deltafactor  df_prime  gbrv_bcc  gbrv_fcc  symbol   md5
            H-high.psp8  0.074830     1.258340  0.028904  0.024726  H        5863396c90149cbe12af496141bde0d0
            ...

        where `gbrv_bcc` and `gbrv_fcc` are the relative errors (in percentage) wrt the AE calculations.
        """
        _TRIALS2KEY = {
            "deltafactor": "dfact_meV",
            "gbrv_bcc": "a0_rel_err",
            "gbrv_fcc": "a0_rel_err",
        }

        rows = []
        for p in self:
            # Extract the dojo_report
            report = p.dojo_report

            row = dict(basename=p.basename, symbol=p.symbol, md5=p.md5)

            for trial, key in _TRIALS2KEY.items():
                # Get results as function of ecut
                try:
                    data = report[trial]
                except KeyError:
                    print("%s does not have %s" % (p.basename, trial))
                    continue

                # Extract the value with highest ecut.
                high_ecut = list(data.keys())[-1]
                row.update({trial: data[high_ecut][key]})
                if trial == "deltafactor":
                    row.update(dict(df_prime=data[high_ecut]["dfactprime_meV"]))

            rows.append(row)

        from pandas import DataFrame
        return DataFrame(rows)

    @add_fig_kwargs
    def plot_dfgbrv_dist(self, **kwargs):
        """
        Plot four distribution plots for the deltafactor, deltafactor prime and the
        relative errors for the GBRV fcc, bcc structures.

        Return: `matplotlib` figure.
        """
        import matplotlib.pyplot as plt
        fig, ax_list = plt.subplots(nrows=2, ncols=2, squeeze=True)
        ax_list = ax_list.ravel()

        frame = self.get_dfgbrv_dataframe()

        import seaborn as sns
        for ax, col in zip(ax_list.ravel(), ["deltafactor", "gbrv_fcc", "df_prime", "gbrv_bcc"]):
            values = frame[col].dropna()
            #print(type(values))
            sns.distplot(values, ax=ax, rug=True, hist=True, kde=False, label=col, bins=kwargs.pop("bins", 50))

            # Add text with Mean or (MARE/RMSRE)
            text = []; app = text.append
            if col in ("deltafactor", "df_prime"):
                app("Mean = %.2f" % values.mean())
            else:
                app("MARE = %.2f" % values.abs().mean())
                app("RMSRE = %.2f" % np.sqrt((values**2).mean()))

            ax.text(0.8, 0.8, "\n".join(text), transform=ax.transAxes)

        return fig
