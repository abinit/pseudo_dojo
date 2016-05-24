"""
Functions providing access to file data. Mainly used to build the public APIs and write unit tests.
"""
from __future__ import print_function, division, unicode_literals

import os

here = os.path.dirname(__file__)

DOJOTABLE_BASEDIRS = [
    "ONCVPSP-PBE-DEV",
    "ONCVPSP-PBE-PDv0.2",
    "ONCVPSP-PBE-PDv0.3",
    "ONCVPSP-PW-DEV"
]


def dojotable_absdir(basedir):
    """
    Return the absolute dirpath of the table from its basename
    """
    if basedir not in DOJOTABLE_BASEDIRS:
        raise RuntimeError(
           "%s is not registered in DOJOTABLE_BASEDIRS\n" 
           "Change pseudo_dojo/pseudos/__init__.py" % basedir) 

    return os.path.join(here, basedir)
