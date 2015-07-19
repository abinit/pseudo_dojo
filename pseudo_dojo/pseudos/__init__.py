"""
Functions providing access to file data.
Mainly used to build the public APIs and write unit tests.
"""
from __future__ import print_function, division, unicode_literals

import os

here = os.path.dirname(__file__)

DOJO_BASEDIRS = [
    "ONCVPSP-PBE",
    "ONCVPSP-PBE-MC2",
]


def dojo_absdir(dojo_basedir):
    """
    Return the absolute path of the directory with the pseudos from its basename
    """
    assert dojo_basedir in DOJO_BASEDIRS
    return os.path.join(here, dojo_basedir)


#def find_ncfiles(top):
#    """
#    Find all netcdf files starting from the top-level directory top.
#    Filenames must be unique. Directories whose start with "tmp_" are
#    excluded from the search.
#
#    Returns:
#        dictionary with mapping: basename --> absolute path.
#    """
#    SILENT = 0
#    ncfiles = {}
#    for dirpath, dirnames, filenames in os.walk(top):
#
#        if "tmp_" in dirpath:
#            continue
#
#        for basename in filenames:
#            apath = os.path.join(dirpath, basename)
#            if basename.endswith(".nc"):
#
#                if basename in ncfiles:
#                    err_msg =  "Found duplicated basename %s\n" % basename
#                    err_msg += "Stored: %s, new %s\n" % (ncfiles[basename], apath)
#
#                    if not SILENT:
#                        import warnings
#                        warnings.warn(err_msg)
#                        #raise ValueError(err_msg)
#                        SILENT += 1
#
#                else:
#                    ncfiles[basename] = apath 
#
#    return ncfiles
