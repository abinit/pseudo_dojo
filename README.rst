
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

    >>> import pseudo_dojo as dojo

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

Global long-term objectives:
----------------------------

    #. To have, on the Web, sets of validated pseudopotentials, for the whole periodic table,
       for different exchange-correlation functionals, with different possibilities of 
       semi-core electrons (e.g. for GW), different cut-off radii (e.g. high pressure application), 
       with an optimal cut-off energy.

    #. To have a Web portal for the generation/validation of new pseudopotentials.

Intermediate objectives (more and more difficult steps):
--------------------------------------------------------

    #. Robust generation of one pseudopotential, be given the atomic number
       (for a given exchange-correlation functional, definition of semi-core electrons, 
       definition of cut-off radii). The cut-off energy does not matter at this stage.

    #. Computation of selected physical properties for selected systems, associated with one 
       given pseudopotential (automatic computation of the cut-off energy, computation of the 
       total energy, the interatomic distance, the lattice parameter of the elemental solid 
       and one oxide, also dimer). Results presented on the Web.

    #. Validation of pseudopotentials with respect to a reference.

    #. Semi-automatic improvement of the generation of the pseudopotential:
       more accurate physical properties, better energy cut-off.

    #. Automatic generation of one set of pseudopotentials, and associated automatic procedure 
       of calculation, validation.

    #. The transfer to the Web of the different results of the objectives.

Components (database):
----------------------

    #. A set of repositories (one for each atomic number). Repository is placed in the 
       directory <Z-symbol> where Z is the atomic  number (three digits) followed by the symbol, 
       with a dash in between, e.g. 001-H, or 092-U.
       Each of these repository is an "ATOM repository", and will contain subdirectories of two 
       kinds, see sections below.

    #. Structuration inside one ATOM repository:
 
         - A ReferenceData subdirectory, that is not tight to one pseudopotential.

         - For each pseudopotential:  
          /[PAW|NC|...]/num_valence_electrons/xc_type/ID , which might called a "Pseudo-atom Box" or patbox.
    
       The num_valence_electrons might be 4e, or 22e ...
       The xc_type might be GGA-PBE or LDA-PW91 or the libcx ID: ex. LIBXC-X001-C012.
       ID will be a digit, e.g. 1, or 2, etc ...
       These IDs will not have any predefined meaning. Some of the pseudo-atom boxes might be good
       for a specific purpose (e.g. GW or High-pressure), but this will be determined from the
       database of results for each pseudo-atom box, by a script, at the demand of one user 
       (or for populating a Web page).

    #. Content of the ReferenceData subdirectory of the "ATOM repository":

        - A (xml, json?) file with the atomic configuration for each possible 
          num_valence_electrons, and other data needed for pseudo-atom generators that are not 
          specific to a pseudo-atom generator. Standard name: atomic_config.xml.

        - Possibly, some CIF files for a elemental solid (or more than one), and for oxide(s) 
          or hydride(s), or potassium-based compounds.

        - A set of master data file (xml) containing the description of the different test systems 
          for the specific atom. Each test system belongs to a test system class: 

             - atom

             - dimer

             - elemental 

             - oxide 

             - hydride or potassium-based compound.

          Within each class it is labelled with an ID (number starting from 1). 
          This description contains insulator/metal and magnetism information, and either
          the name of the cif file to be used, or the reference length for the dimer.
          Standard name: <Z-symbol>.description_<class>_ID.xml

    #. Content of one pseudo-atom Box:

       - Subdirectories: 

                - <name_of_generator>

                - atom_X, dimer_X

                - elemental_X, oxide_X

                - possibly hydride_X or K_X . 

         Where <name_of_generator> might be atompaw, or ape, or fhi98pp, or ...
         And where X is the ID defined in D3c.

        - A (xml) summary file containing metadata concerning this pseudo-atom box, obtained by 
          running the applications in the different subdirectories, and also describing the 
          validation criteria (this implies a set of runs).
          Standard name: <Z-symbol>.summary.[PAW|NC|...].#valence_electrons.<xc_type>.ID.xml

    #. Content of the <name_of_generator> subdirectory of the pseudo-atom Box:

       - Optionally, the specific input data needed for the generator (PAW or NC), to complement 
         the content of the atomic_config file. Typically cut-off radii.
         Standard name atomic_data_<name_of_generator>.xml, e.g. atomic_data_atompaw.xml
      
       - A pseudo-atom data generator input file (PAW or NC) - might have been automatically generated
         from atomc_config and the file in D5a. Standard name <name_of_generator>.in

       - A pseudo-atom data file (PAW or NC) - has been automatically generated (output of the atomic generator).
         Standard name <Z-symbol>.pseudoatom.[PAW|NC|...].#valence_electrons.<xc_type>.ID. 

       - <standard_postfix_for_the_generator>
         This is the pseudopotential file, or the PAW atomic data file.
         The <standard_postfix_for_the_generator> might be .fhi or .pawps , or other postfix.

    #. Content of the <class>_X (where class is atom, dimer, elemental, oxide, hydride ...)
       D6a Subdirectories abinit_runY and elk_runY, where Y is an integer starting from 1.

    #. Content of the abinit_runY directory
       D7a This is a working directory for one abinit run. It contains an ABINIT input file 
       usually automatically generated from D3c, specialized for the pseudo-atom box and the system.
       D7b For Y=1 : determination of a basic k point grid, using kptrlen and prtkpt. 
       Can be used by elk, see elk_1.
       D7c For Y=2 : computation of total energy as a function of ecut, for basic k point grid, 
       and, for metals, using the tsmear determined by elk_1.

    #. Content of the elk_runY subdirectory:

       - This is a working directory for one elk run. It contains an ELK input file 
         usually automatically generated from D3c, specialized for the pseudo-atom box 
         and the system, and using the k point grid determined by D7b .

       - For Y=1, determination of the tsmear.

Components (software).
----------------------

They should be placed inside ABINIT package psps/script, for testing/coherency purposes 
across the different <Z-symbol> directories.

   #. A "pseudo-atom box" creator (init_patbox.py), to be called inside the psps/<Z-symbol> directory.
      (propose options for the path described in D2, then create the path, 
      and the directories of the pseudo-atom box, and also bzr add the dirs)

   #. A cif2cml translator, to go from D3b to D3c.

   #. A script to initialize the file <name_of_generator>.in mentioned in D5b from D3a 
      atomic_config.xml and D5a atomic_data_<name_of_generator>.xml
 
   #. A pseudopotential generator, e.g. ATOMPAW  (already placed inside the ABINIT package)

   #. A driver of abinit: generation of abinit input files, running of abinit, gathering of the data in D4b. 
      ACTUALLY NEED A LIST OF TASKS / VALIDATION CRITERIA /  to be defined.

   #. A driver of elk, and a binary for elk.
   
   #. A validator.

   #. More scripts to be added ...

Miscellaneous
-------------

  #. Reference oxygen PAW data files for different XC functionals, reference hydrogen PAW data files
     for different XC functionals.  Placed in the abinit/psps/RefPseudoAtoms subdirectory of 
     the ABINIT package. And to be copied in the patbox at init time.

Strategy
--------

   - Work component by component, by placing these components under version management and 
     automatic testing, with appropriate hardware.

   - Define the files and their format (including metadata) in an iterative way, with possibilities 
     to regenerate them in an automatic way

   - Gradual understanding of the CPU constraints, memory constraint, and human time needed.

   - Adjust the objectives to stay realist.

TO BE KEPT IN MIND FOR FURTHER SPECIFICATION
--------------------------------------------

   - set up of a bot (on the machine nazgul): be given the ABINIT branch, and the pseudopotentials =>
     computation of the physical characteristics of this pseudopotential

   - set up of the corresponding "on-demand" mechanism

   - set up a new waterfall: the list of files that will be provided will be quite different 
     from the usual bots

   - set up of a new Web window to visualize the files (to be discussed).


License
=======

pseudo_dojo is released under the GPL License. The terms of the license are as follows:

pseudo_dojo is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

pseudo_dojo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with pseudo_dojo.  
If not, see <http://www.gnu.org/licenses/>.

