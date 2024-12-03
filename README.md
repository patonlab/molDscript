![molDscript](https://github.com/patonlab/molDscript/blob/main/moldscript/molDscript.png)
===

[![CircleCI](https://dl.circleci.com/status-badge/img/circleci/JDvVi58JeRw4LYzfeJGsjn/T7u88FmqfkcE7c6vKveH1L/tree/main.svg?style=shield&circle-token=CCIPRJ_3nGjXb4n3dHaAo6mQ67TBk_5ce95f5de89641ed836cbe55488e9b11f28c43d3)](https://dl.circleci.com/status-badge/redirect/circleci/JDvVi58JeRw4LYzfeJGsjn/T7u88FmqfkcE7c6vKveH1L/tree/main)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


A Python package for ...

## Requirements

* cclib
* DBSTEP
* rdkit
* yaml

## Key parameters

* --varfile (str): Specify a .txt file of arguments with the following format:
```
parameter1: value1
parameter2: value2
``` 
* --opt (str): Specify the path to a folder containing your optimization files
* --nbo (str): Specify the path to a folder containing your nbo files
* --nmr (str): Specify the path to a folder containing your nmr files
* --fukui_neutral (str): Specify the path to a folder containing your neutral fukui files
* --fukui_reduced (str): Specify the path to a folder containing your reduced fukui files
* --fukui_oxidized (str): Specify the path to a folder containing your oxidized fukui files
* --charges (str): Specify the path to a folder containing your charges files
* --fmo (str): Specify the path to a folder containing your fmo and dipole files
* --ad_reduced (str): Specify the path to a folder containing your adiabatic reduced files
* --ad_oxidized (str): Specify the path to a folder containing your adiabatic oxidized files
* --substructure ("str"): specify the substructure you want to search for in the molecule
* --volume : Indicate you want the buried volume of the atoms in the substructure match
* --vall : Indicate you want the volume of all atoms in every molecule
* --radius (float or [float, float]): Speficy the radius/radii you want buried volume to be calculated for
* --boltz : Indicate you want the boltzmann weighted average of the conformers
* --temp (float): Indicate the temperature to use for boltzmann weighting
* --min_max : Indicate you want the minimum, maximum, and range values of the parameters

# Installation

## pip

```shell
$ pip install XXX
```

## conda

```shell
$ conda install -c conda-forge XXX
```

# Usage

XXX can be imported as a Python module that is easily integrated into
workflows. Here is an example for ...

```shell
>>> from XXX import YYY
>>> etc
```

It can also be used from the command line.

```console
$ python -m moldscript --varfile arguments.txt
```
For further information, see the separate [documentation](https:/XXX).

# Packages Supported

* Gaussian
* Orca
* xTB

# Features

### Molecule level parameters

* Dipole and quadrupole moments
* Molecular apolarizabilities and hyperpolarizabilities
* Orbital energies: HOMO energy, LUMO energy, HOMO-LUMO gap
* Global reactivity descriptors: electronegativity, hardness, softness, global electrophilicity index
* Ionization Potential and Electron Affinity
* Molecular volume, Van der Waals volume, Solvent-accessible volume (not yet implemented)
* Surface Area: Van der Waals surface area, Solvent-accessible surface area (SASA), Polar surface area (PSA) (not yet implemented)

### Bond level parameters

* Bond Orders: Wiberg bond indices

### Atom level parameters

* Partial charges: Natural Population Analysis (NPA) charges, Hirshfeld charges
* Atomic apolarizabilities and hyperpolarizabilities (not yet implemented)
* Fukui indices: Nucleophilic fukui function (f-), Electrophilic fukui function (f+), Radical fukui function (f0)
* NMR Chemical Shifts
* Atomic Volumes (not yet implemented)
* Buried Volumes (not yet implemented)
* Vol2Vec parameters (not yet implemented)


### Testing
details 

Tests can be run with the `pytest -v` command. There are a number of additional command line arguments to explore.

### Documentation
something about readthedocs


# Acknowledgements

This work was carried out in the [Paton Laboratory at Colorado State University](https://patonlab.colostate.edu), supported by the [NSF Center for Computer-Assisted Synthesis](https://ccas.nd.edu/), grant number [CHE-1925607](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2202693&HistoricalAwards=false)

In particular, the following people have contributed significantly to developing its functionality:

* [Shree Sowndarya](https://github.com/shreesowndarya)
* [Jake King](https://github.com/j77king)
* [Robert Paton](https://github.com/bobbypaton)

