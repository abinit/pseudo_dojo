"""
This module provides a database with the acell lattice parameter for Rocksalt Nitrides
wit Rare-earth cations
Client code should use the official API df_database(xc) to access the database.

Example::

    df = raren_database("PBE")
    df
"""
from __future__ import print_function, division, unicode_literals

import os
import pandas as pd

from pseudo_dojo.core.pseudos import Pseudo

_ROOT = os.path.dirname(__file__)


def _get_database(xc):
    # data_RE_nitrides_PBE.txt
    path = os.path.join(_ROOT, "data", 'data_RE_nitrides_%s.txt' % str(xc))
    if not os.path.isfile(path):
        raise ValueError("No such file: %s" % path)

    with open(path, 'rt') as f:
        lines = f.readlines()

    # acell vs RE cations Rocksalt Nitrides
    # atom Z  exp   VASP   ONCVPSP VLab    JTH  FOO
    #
    #La 57 5.29 5.265 5.252 na na  56
    data = []
    for i, line in enumerate(lines):
        line = line.strip().lstrip()
        if line.startswith("#"): line = line[1:]
        if i == 0:
            continue
        elif i == 1:
            columns = [t for t in line.split(" ") if t]
        else:
            data.append([t for t in line.split(" ") if t])

    #print(len(data), len(columns))
    table = pd.DataFrame(data, columns=columns)
    table.set_index('atom', inplace=True)
    table.index.name = None
    table = table.apply(pd.to_numeric, errors='coerce')

    # Compute reference value here.
    #table['ref'] = table[['VASP', 'VLab']].mean(1)
    table['ref'] = table['VASP']
    return RareEarthNDatabase(xc, table)


class RareEarthNDatabase(object):
    """
    This object stores the results of the structural relaxation
    """
    def __init__(self, xc, table):
        self.xc, self.table = xc, table

    def get_n_pseudo(self, pseudo):
        if pseudo.ispaw:
            raise NotImplementedError()
        else:
            # N_PBE.psp8
            n_pseudo = os.path.join(_ROOT, "data", "N_%s.psp8" % str(self.xc))

        n_pseudo = Pseudo.from_file(n_pseudo)
        assert n_pseudo.xc == pseudo.xc
        return n_pseudo

##########################################################################################
# Official API to access the database.
##########################################################################################


# Mapping XC --> Database.
__RAREN_DATABASE_XC = None


def raren_database(xc):
    """
    Returns the database with the reference results associated to this XC functional
    xc is the exchange-correlation functional e.g. PBE, PW
    """
    global __RAREN_DATABASE_XC
    if __RAREN_DATABASE_XC is None:
        __RAREN_DATABASE_XC = {}

    # Create xc database and cache it.
    if xc not in __RAREN_DATABASE_XC:
        db = _get_database(xc=xc)
        __RAREN_DATABASE_XC[xc] = db

    return __RAREN_DATABASE_XC[xc]
