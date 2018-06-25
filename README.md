
PseudoDojo is an open-source Python framework for generating and validating pseudopotentials.
It uses the [AbiPy](https://github.com/abinit/abipy) package, for developing and systematically 
testing pseudopotentials. 
At present it contains seven different batteries of tests executed with [ABINIT](https://github.com/abinit/abinit) 
([delta-gauge](https://molmod.ugent.be/deltacodesdft), 
revised delta gauge by [Jollet et al](http://www.sciencedirect.com/science/article/pii/S0010465513004359), 
[GBRV](https://www.physics.rutgers.edu/gbrv/) tests for fcc/bcc/compounds, 
phonons at the Gamma point and tests for the presence of ghost-states below and above the Fermi level).

The tests are performed as a function of the energy cutoff (with the exception of ghosts and compounds) 
and the results are then used to provide hints for the energy cutoff for actual production calculations. 
We keep track of the results of the various validation tests for each pseudopotential with the `djrepo` file, 
a text document in JSON format produced by the python code at the end of the test. 
For further details, please consult our recent [cond-mat paper](https://arxiv.org/abs/1710.10138).

The PseudoDojo code is hosted on [github](https://github.com/abinit/pseudo_dojo) 
but we also provide a user web-interface at <http://www.pseudo-dojo.org>
There you will find pseudopotential files that can be used immediately, 
as well as the corresponding input files if you need to change or tune or change some parameters 
(e.g. the XC functional).

The pseudopotential files are available on the web-site in the ABINIT `psp8` format, 
in the `UPF2` format and in the `PSML` 1.1 XML format shared by SIESTA and ABINIT. 
The input files, the results of the generation, and the test results are presented via jupyter notebooks 
as static HTML pages. 
One can hence easily compare pseudopotentials for a given element and then select the most appropriate 
one according to a chosen criterion (e.g. efficiency vs accuracy).  

How to cite the PseudoDojo project
----------------------------------

If you use the PseudoDojo pseudopotentials in your research, 
please consider citing the works mentioned in [this page](http://www.pseudo-dojo.org/faq.html).

Getting PseudoDojo
------------------

Developmental version
---------------------

The developmental version is at the PseudoDojo's `Github repo <https://github.com/gmatteo/pseudo_dojo>`_. 
After cloning the source, you can type::

    python setup.py install

or, alternatively,

    python setup.py develop

to install the package in developmental mode.

<!--
Stable version
==============

The version at the Python Package Index (PyPI) is always the latest stable
release that will be hopefully, be relatively bug-free. The easiest way to
install PseudoDojo is to use pip, as follows::

    pip install pseudo_dojo
-->


Project PseudoDojo: global view
===============================

Global long-term objectives:
----------------------------

- To have, on the Web, sets of validated pseudopotentials, for the whole periodic table,
  for different exchange-correlation functionals, with different possibilities of 
  semi-core electrons (e.g. for GW), different cut-off radii (e.g. high pressure application), 
  with an optimal cut-off energy.

- To have a Web portal for the generation/validation of new pseudopotentials.

- Computation of selected physical properties for selected systems, associated with one 
  given pseudopotential (automatic computation of the cut-off energy, computation of the 
  total energy, the interatomic distance, the lattice parameter of the elemental solid 
  and one oxide, also dimer). Results presented on the Web.

- Validation of pseudopotentials with respect to a reference.

License
=======

The pseudopotential files as well as the input files are released under the 
[CC BY 4.0 license](https://creativecommons.org/licenses/by/4.0/legalcode)
This license lets you distribute, remix, tweak, and build upon our work, even commercially, 
as long as you credit the PseudoDojo project for the original creation. 

The PseudoDojo code (the python code used to generate/validate the pseudos) 
is released under the GPL License. 
The terms of the license are as follows:

PseudoDojo is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

PseudoDojo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PseudoDojo.  
If not, see <http://www.gnu.org/licenses/>.
