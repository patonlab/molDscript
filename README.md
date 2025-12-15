![molDscript](https://github.com/patonlab/molDscript/blob/main/moldscript/molDscript.png)
===

[![CircleCI](https://dl.circleci.com/status-badge/img/circleci/JDvVi58JeRw4LYzfeJGsjn/T7u88FmqfkcE7c6vKveH1L/tree/main.svg?style=shield&circle-token=CCIPRJ_3nGjXb4n3dHaAo6mQ67TBk_5ce95f5de89641ed836cbe55488e9b11f28c43d3)](https://dl.circleci.com/status-badge/redirect/circleci/JDvVi58JeRw4LYzfeJGsjn/T7u88FmqfkcE7c6vKveH1L/tree/main)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

MolDscript is a Python workflow that converts Density Functional Theory (DFT) and related quantum chemistry outputs into descriptor tables ready for machine learning or benchmarking. It wraps `cclib`, `RDKit`, `DBSTEP`, and `pandas` to align Gaussian, ORCA, and xTB calculations and write consistent molecule-, bond-, and atom-level CSV files.

## Highlights
- Parse optimization, single-point, NBO, NMR, charge, FMO, and Fukui calculations without manual file editing.
- Match conformer ensembles, apply SMARTS-based substructure filters, and compute DBSTEP buried volumes on demand.
- Generate ensembles (Boltzmann weighted, min/mnax within population windows, lowest-energy snapshots) in a single run.
- Emit descriptor CSVs alongside module logs (`MOLDSCRIPT_*.dat`) for traceability.

## Installation
```shell
git clone https://github.com/patonlab/molDscript.git
cd molDscript
pip install -e .
```
Open Babel (optional) can be installed from conda-forge:
```shell
conda install -c conda-forge openbabel
```

## Quick Start
```shell
python -m moldscript \
  --opt calculations/opt --fmo calculations/fmo 
  --suffix_fmo fmo_suffix --suffix_nbo nbo_suffix
  --nbo calculations/nbo  

```
Prefer storing options in a key:value text file? Use `--varfile inputs.txt`; command-line flags override values loaded from the file.

## Core Inputs & Flags
- `--opt PATH` (required) - baseline optimization files and conformer metadata.
- `--spc PATH` - single-point energies that replace optimization SCF energies.
- `--nbo`, `--nmr`, `--charges`, `--fmo` PATH - add module-specific descriptors; pair with `--suffix_*` when filenames include extra tokens.
- `--fukui_neutral`, `--fukui_reduced`, `--fukui_oxidized` PATH - supply all three charge states for vertical IE/EA and condensed Fukui functions.
- `--substructure SMARTS` - limit atom/bond descriptors to a SMARTS match; combine with `--volume` or `--vall` and optional `--radius` list for DBSTEP buried volumes.
- `--boltz`, `--min_max`, `--lowe` - compute Boltzmann-weighted averages, min/max/range tables (using `--cut`), and lowest-energy snapshots. Adjust `--temp` (K) as needed.
- `--output PREFIX` - prepend every generated filename; append a slash to target a directory. Use `--no_mol`, `--no_atom`, `--no_bond`, or `--no_bond_filter` to tailor CSV output.

## Output Artefacts
- `molecule_level.csv`, `bond_level.csv`, `atom_level.csv` - aligned descriptors per calculation, bond pair, or atom.
- `ensemble_*.csv`, `boltzmann_weights.csv` - created when `--boltz` is enabled.
- `min_max_range_*.csv`, `lowest_energy_*.csv` - created when `--min_max` or `--lowe` are requested.
- `MOLDSCRIPT_*.dat` - per-module logs capturing provenance and CPU-time summaries.

## Documentation
The Read the Docs site (coming soon) will provide the full user guide: [https://moldscript.readthedocs.io](https://moldscript.readthedocs.io)

## Dependencies
Key Python dependencies include `pandas`, `cclib` (latest GitHub version for the most up-to-date package compatability), `dbstep`, `rdkit`, `networkx`, `numpy`, and `periodictable`.

## Supported Quantum Packages
- Gaussian
- ORCA
- xTB (optimizations)

## Testing
Run `pytest -v` from the project root to execute the test suite.

## Acknowledgements
This work was carried out in the [Paton Laboratory at Colorado State University](https://patonlab.colostate.edu), supported by the [NSF Center for Computer-Assisted Synthesis](https://ccas.nd.edu/) (grant [CHE-1925607](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2202693&HistoricalAwards=false)).

Contributors include [Shree Sowndarya](https://github.com/shreesowndarya), [Jake King](https://github.com/j77king), and [Robert Paton](https://github.com/bobbypaton).
