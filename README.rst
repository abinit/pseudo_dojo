
Pseudo_dojo is a open-source Python framework for generaring and validating  pseudo potentials.

pseudo_dojo is free to use. However, we also welcome your help to improve this
library by making your own contributions.  These contributions can be in the
form of additional tools or modules you develop, or even simple things such
as bug reports. Please report any bugs and issues at pseudo_dojo's `Github page
<https://github.com/gmatteo/pseudo_dojo>`_. 

Getting pseudo_dojo
================

Stable version
--------------

The version at the Python Package Index (PyPI) is always the latest stable
release that will be hopefully, be relatively bug-free. The easiest way to
install pseudo_dojo on any system is to use easy_install or pip, as follows::

    easy_install pseudo_dojo

or::

    pip install pseudo_dojo

Developmental version
---------------------

The developmental version is at the pseudo_dojo's `Github repo <https://github.com/gmatteo/pseudo_dojo>`_. 
After cloning the source, you can type::

    python setup.py install

or to install the package in developmental mode::

    python setup.py develop

Requirements
============

All required dependencies should be automatically taken care of if you
install pseudo_dojo using easy_install or pip. Otherwise, these packages should
be available on `PyPI <http://pypi.python.org>`_.

1. Python 2.7+ required. 
2. numpy 
3. scipy 0.10+ 
4. matplotlib 1.1+ 
5. periodictable

Using pseudo_dojo
=================

Basic usage
-----------

Useful aliases for commonly used objects are now provided, similar in style to
numpy. Supported objects include Element, Composition, Structure, Molecule,
Spin and Orbital. Here are some quick examples of the core capabilities and
objects:

.. code-block:: pycon

    >>> import pseudo_dojo as mg
    >>>
    >>> si = mg.Element("Si")
    >>> si.atomic_mass
    28.0855
    >>> si.melting_point
    u'1687 K'
    >>>
    >>> comp = mg.Composition("Fe2O3")
    >>> comp.weight
    159.6882
    >>> #Note that Composition conveniently allows strings to be treated just
    >>> #like an Element object.
    >>> comp["Fe"]
    2.0
    >>> comp.get_atomic_fraction("Fe")
    0.4
    >>> lattice = mg.Lattice.cubic(4.2)
    >>> structure = mg.Structure(lattice, ["Cs", "Cl"],
    ...                          [[0, 0, 0], [0.5, 0.5, 0.5]])
    >>> structure.volume
    74.088000000000008
    >>> structure[0]
    PeriodicSite: Cs (0.0000, 0.0000, 0.0000) [0.0000, 0.0000, 0.0000]
    >>>
    >>> #Integrated symmetry tools from spglib.
    >>> from pseudo_dojo.symmetry.finder import SymmetryFinder
    >>> finder = SymmetryFinder(structure)
    >>> finder.get_spacegroup_symbol()
    'Pm-3m'
    >>>
    >>> #Writing out a POSCAR file for VASP calculations.
    >>> poscar = Poscar(structure)
    >>> mg.write_structure(structure, "POSCAR")
    >>>
    >>> #Reading a structure from a file.
    >>> structure = mg.read_structure("POSCAR")

The above illustrates only the most basic capabilities of pseudo_dojo.

.. note:: Examples

    A good way to explore the functionality of pseudo_dojo is to look at examples.
    We have created a `Github wiki page
    <https://github.com/materialsproject/pseudo_dojo/wiki>`_ to allow users to
    share their Github gists (essentially mini git repos of scripts)
    performing various kinds of functions with pseudo_dojo. Please feel free to
    check them out and we welcome your contributions as well!

===============================
Project PseudoDojo: global view
===============================

.. include:: PROJECT.rst

License
=======

pseudo_dojo is released under the GPL License. The terms of the license are as follows:

.. literalinclude:: LICENSE

