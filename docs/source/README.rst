
.. image:: ../../moldscript/molDscript.png
    :alt: molDscript

.. image:: https://dl.circleci.com/status-badge/img/circleci/JDvVi58JeRw4LYzfeJGsjn/T7u88FmqfkcE7c6vKveH1L/tree/main.svg?style=shield&circle-token=CCIPRJ_3nGjXb4n3dHaAo6mQ67TBk_5ce95f5de89641ed836cbe55488e9b11f28c43d3
    :target: https://dl.circleci.com/status-badge/redirect/circleci/JDvVi58JeRw4LYzfeJGsjn/T7u88FmqfkcE7c6vKveH1L/tree/main

.. image:: https://img.shields.io/badge/License-MIT-yellow.svg
    :target: https://opensource.org/licenses/MIT


|
Introduction
------------
molDscript is a Python package that curates a variety of molecular, bond, and atom-level properties from Density Functional Theory (DFT) output files.
These properties are curated into .csv files that can be easily used for machine learning, benchmarking, and other computational chemistry applications.

Installation
------------

First, clone the repository from https://github.com/patonlab/molDscript.git. Then, either manually add the path to your PYTHONPATH or install the package using pip while in the cloned repository.

.. code-block:: shell

     $ pip install .

Requirements
------------

* cclib
* DBSTEP
* rdkit
* yaml

Using molDscript
--------------

.. code-block:: text

     python -m moldscript [--varfile /path/to/arguments.txt] [--opt /path/to/optimization/files] [--nbo /path/to/nbo/files] 
     [--nmr /path/to/nmr/files] [--fukui_neutral /path/to/neutral_fukui/files] [--fukui_reduced /path/to/reduced_fukui/files] 
     [--fukui_oxidized /path/to/oxidized_fukui/files] [--charges /path/to/charges/files] [--fmo /path/to/fmo/files] 
     [--ad_reduced /path/to/reduced_adiabatic/files] [--ad_oxidized /path/to/oxidized_adiabatic/files] [--substructure "substructure"] 
     [--volume] [--vall] [--radius 1.0] [--boltz] [--temp 298.15] [--min_max] 

* ``--varfile``  Specify a .txt file of arguments with the following format::

.. code-block:: text

     parameter1: value1
     parameter2: value2

* ``--opt`` Specify the path to a folder containing your optimization files. This is required for the program to run
* ``--nbo`` Specify the path to a folder containing your nbo files
* ``--nmr`` Specify the path to a folder containing your nmr files
* ``--fukui_neutral`` Specify the path to a folder containing your neutral fukui files
* ``--fukui_reduced`` Specify the path to a folder containing your reduced fukui files
* ``--fukui_oxidized`` Specify the path to a folder containing your oxidized fukui files
* ``--charges`` Specify the path to a folder containing your charges files. This will default to the spc files (if available) or opt files if not specified.
* ``--fmo`` Specify the path to a folder containing your fmo and dipole files. This will default to the spc files (if available) or opt files if not specified.
* ``--ad_reduced`` Specify the path to a folder containing your adiabatic reduced files
* ``--ad_oxidized`` Specify the path to a folder containing your adiabatic oxidized files
* ``--substructure`` Specify the substructure you want to search for in the molecule
* ``--volume`` Indicate you want the buried volume of the atoms in the substructure match
* ``--vall`` Indicate you want the volume of all atoms in every molecule
* ``--radius`` Specify the radius/radii you want buried volume to be calculated for
* ``--boltz`` Indicate you want the boltzmann weighted average of the conformers
* ``--temp`` Indicate the temperature to use for boltzmann weighting
* ``--min_max`` Indicate you want the minimum, maximum, and range values of the parameters


Packages Supported
------------------

* Gaussian
* Orca
* xTB (for opt files only)

Features
--------

Molecule level parameters
~~~~~~~~~~~~~~~~~~~~~~~~~

* Dipole and quadrupole moments
* Molecular polarizabilities and hyperpolarizabilities
* Orbital energies: HOMO energy, LUMO energy, HOMO-LUMO gap
* Global reactivity descriptors: electronegativity, hardness, softness, global electrophilicity index
* Ionization Potential and Electron Affinity
* Molecular volume, Van der Waals volume, Solvent-accessible volume (not yet implemented)
* Surface Area: Van der Waals surface area, Solvent-accessible surface area (SASA), Polar surface area (PSA) (not yet implemented)

Bond level parameters
~~~~~~~~~~~~~~~~~~~~~

* Bond Orders: Wiberg bond indices

Atom level parameters
~~~~~~~~~~~~~~~~~~~~~

* Partial charges: Natural Population Analysis (NPA) charges, Hirshfeld charges
* Atomic polarizabilities and hyperpolarizabilities (not yet implemented)
* Fukui indices: Nucleophilic fukui function (f-), Electrophilic fukui function (f+), Radical fukui function (f0)
* NMR Chemical Shifts
* Atomic Volumes (not yet implemented)
* Buried Volumes (not yet implemented)
* Vol2Vec parameters (not yet implemented)

Testing
-------

Tests can be run with the `pytest -v` command. 



Acknowledgements
----------------

This work was carried out in the `Paton Laboratory at Colorado State University <https://patonlab.colostate.edu>`_, supported by the `NSF Center for Computer-Assisted Synthesis <https://ccas.nd.edu>`_, grant number `CHE-1925607 <https://www.nsf.gov/awardsearch/showAward?AWD_ID=2202693&HistoricalAwards=false>`_.

In particular, the following people have contributed significantly to developing its functionality:

* `Shree Sowndarya <https://github.com/shreesowndarya>`_
* `Jake King <https://github.com/j77king>`_
* `Robert Paton <https://github.com/bobbypaton>`_
