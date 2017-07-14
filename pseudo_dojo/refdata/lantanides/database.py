"""
This module provides a database with the acell lattice parameter for Rocksalt Nitrides
wit Rare-earth cations
Client code should use the official API df_database(xc) to access the database.

Example::

    df = renrs_database("PBE")
    df
"""
from __future__ import print_function, division, unicode_literals

import os
import pandas as pd


def _get_dataframe(xc):
    path = os.path.join(os.path.dirname(__file__), 'data_RE_nitrides_%s.txt' % str(xc))
    with open(path, 'rt') as f:
        lines = f.readlines()

    # acell vs RE cations Rocksalt Nitrides
    # atom Z  exp   VASP   ONCVPSP VLab    JTH  FOO
    #
    #La 57 5.29 5.265 5.252 na na  56
    data = []
    for i, line in enumerate(lines):
        if i == 1:
            columns = []
            for item in line.split(' '):
                if item not in ['', '#', '\n']:
                    columns.append(item)
        if '#' not in line:
            data.append(line.strip().split(' '))

    df = pandas.DataFrame(data, columns=columns)
    df.set_index('atom', inplace=True)
    df.index.name = None
    df = df.apply(pd.to_numeric, errors='coerce')

    df['ref'] = df[['VASP','VLab']].mean(1)
    return df


##########################################################################################
# Official API to access the database.
##########################################################################################

# Mapping XC --> Database.
__REN_DATABASE_XC = None


def renrs_database(xc):
    """
    Returns the database with the reference results associated to this XC functional
    xc is the exchange-correlation functional e.g. PBE, PW
    """
    global __REN_DATABASE_XC
    if __REN_DATABASE_XC is None:
        __REN_DATABASE_XC = {}

    # Create xc database and cache it.
    if xc not in __REN_DATABASE_XC:
        df = _get_dataframe(xc=xc)
        __REN_DATABASE_XC[db.xc] = df

    return __REN_DATABASE_XC[xc]
