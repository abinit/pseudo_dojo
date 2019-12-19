#!/usr/bin/env python
# flake8: noqa
"""Setup script for pseudo_dojo."""
import sys
import os

from glob import glob
from setuptools import find_packages


def file_doesnt_end_with(test, endings):
    """
    A little utility we'll need below, since glob() does NOT allow you to do exclusion on multiple endings!

    Return true if test is a file and its name does NOT end with any
    of the strings listed in endings.
    """
    if not os.path.isfile(test):
        return False
    for e in endings:
        if test.endswith(e):
            return False
    return True

#---------------------------------------------------------------------------
# Basic project information
#---------------------------------------------------------------------------

#with open("README.rst") as f:
#    long_desc = f.read()
#    ind = long_desc.find("\n")
#    long_desc = long_desc[ind + 1:]

# release.py contains version, authors, license, url, keywords, etc.
release_file = os.path.join('pseudo_dojo', 'core', 'release.py')

with open(release_file) as f:
    code = compile(f.read(), release_file, 'exec')
    exec(code)


#def find_packages():
#    """Find all packages."""
#    return find_packages(exclude=())


def find_package_data():
    """Find pseudo_dojo's package_data."""
    # This is not enough for these things to appear in an sdist.
    # We need to muck with the MANIFEST to get this to work
    package_data = {
        'pseudo_dojo.refdata.deltafactor.data': ['*.txt', 'CIFs/*', 'history/*'],
        'pseudo_dojo.refdata.gbrv.data': ['*.csv'],
        'pseudo_dojo.refdata.lantanides.data': ['*'],
        'pseudo_dojo.pseudos': ["ONCVPSP-PBE/*/*",],
    }
    return package_data


def find_scripts():
    """Find pseudo_dojo scripts."""
    scripts = []
    # All python files in abipy/scripts
    pyfiles = glob(os.path.join('pseudo_dojo', 'scripts', "*.py") )
    scripts.extend(pyfiles)

    return scripts


def cleanup():
    """Clean up the junk left around by the build process."""
    if "develop" not in sys.argv:
        import shutil
        try:
            shutil.rmtree('pseudo_dojo.egg-info')
        except:
            try:
                os.unlink('pseudo_dojo.egg-info')
            except:
                pass


# List of external packages we rely on.
# Note setup install will download them from Pypi if they are not available.
with open("requirements.txt", "rt") as fh:
    install_requires = [s.strip() for s in fh]

#---------------------------------------------------------------------------
# Find all the packages, package data, and data_files
#---------------------------------------------------------------------------

# Get the set of packages to be included.
my_packages = find_packages(exclude=())

my_scripts = find_scripts()

my_package_data = find_package_data()

#data_files = find_data_files()

# Create a dict with the basic information
# This dict is eventually passed to setup after additional keys are added.
setup_args = dict(
      name=name,
      version=version,
      description=description,
      long_description=long_description,
      author=author,
      author_email=author_email,
      #url=url,
      #download_url=download_url,
      license=license,
      platforms=platforms,
      keywords=keywords,
      install_requires=install_requires,
      packages=my_packages,
      package_data=my_package_data,
      scripts=my_scripts,
      )


if __name__ == "__main__":
    from setuptools import setup
    setup(**setup_args)
    cleanup()
