# coding: utf-8
"""Release data for the pseudo_dojo project."""

# Name of the package for release purposes.  This is the name which labels
# the tarballs and RPMs made by distutils, so it's best to lowercase it.
name = 'pseudo_dojo'

# version information.  An empty _version_extra corresponds to a full
# release.  'dev' as a _version_extra string means this is a development version
_version_major = 0
_version_minor = 1
_version_micro = '0'  # use '' for first of series, number for 1 and above
#_version_extra = ''  # Uncomment this for full releases

# Construct full version string from these.
_ver = [_version_major, _version_minor]
if _version_micro: _ver.append(_version_micro)
#if _version_extra: _ver.append(_version_extra)

__version__ = '.'.join(map(str, _ver))

version = __version__  # backwards compatibility name

description = "Set of python modules and scripts to automate pseudopotential generation.",

long_description = \
"""
The goal of pseudo_dojo is to create a environment for interactive and exploratory
generation of pseudopotetianls.
To support this goal, pseudo_dojo provides:

* A set of pythons objects that automates the creation of input files.

* A set of scripts for performing common tasks such as plotting the AE, PS wavefunctions.

The latest development version is always available from ...
site <http:// TODO >.
"""

license = 'GPL'

authors = {
          'Matteo': ('Matteo Giantomassi','gmatteo at gmail.com'),
          }

author = 'Matteo Giantomassi'

author_email = 'gmatteo at gmail com'

#url = 'http://abinit.org'
#download_url = 'http://abinit.org/abipy'

platforms = ['Linux', 'darwin']

keywords = ["pseudopotentials", "DFT", "ab-initio", "first principles", "ABINIT"]
