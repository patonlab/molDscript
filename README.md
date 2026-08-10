![molDscript](https://github.com/patonlab/molDscript/blob/main/moldscript/molDscript.png)
===

[![CircleCI](https://dl.circleci.com/status-badge/img/circleci/JDvVi58JeRw4LYzfeJGsjn/T7u88FmqfkcE7c6vKveH1L/tree/main.svg?style=shield&circle-token=CCIPRJ_3nGjXb4n3dHaAo6mQ67TBk_5ce95f5de89641ed836cbe55488e9b11f28c43d3)](https://dl.circleci.com/status-badge/redirect/circleci/JDvVi58JeRw4LYzfeJGsjn/T7u88FmqfkcE7c6vKveH1L/tree/main)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

MolDscript is a Python workflow that converts Density Functional Theory (DFT) and related quantum chemistry outputs into descriptor tables ready for machine learning or benchmarking. It wraps `cclib`, `RDKit`, `DBSTEP`, and `pandas` to align Gaussian, ORCA, and xTB calculations and write consistent molecule-, bond-, and atom-level CSV files.

## Highlights
- Parse optimization, single-point, NBO, NMR, charge, FMO, and Fukui calculations without manual file editing.
- Read CREST and ORCA GOAT multi-conformer XYZ files directly and summarize steric variation without requiring a higher-level calculation.
- Match quantum-chemistry conformer ensembles, apply SMARTS-based substructure filters, and compute DBSTEP buried volumes on demand.
- Generate ensembles (Boltzmann weighted, min/max within population windows, lowest-energy snapshots) in a single run.
- Emit descriptor CSVs alongside a single run audit log (`MOLDSCRIPT.dat`) for traceability.

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

For a geometry-only conformer ensemble:

```shell
python -m moldscript \
  --ensemble conformer_search \
  --ensemble_radii "[3.5]"
```

Directory input recognizes CREST `crest_conformers.xyz` /
`*_crest_conformers.xyz` files and ORCA GOAT `*.finalensemble.xyz` files. Any
equivalent multi-frame XYZ filename can be passed directly. Comment-line
energies are interpreted as Hartree; the minimum is retained internally but is
not added as an ensemble descriptor.

Only compact steric summaries are added to `molecule_level.csv`: min, max, and
range for mass-weighted radius of gyration and relative shape anisotropy.
Buried volume is instead an atom-level descriptor: every atom is used as a
center, and every requested radius adds the minimum, maximum,
Boltzmann-weighted mean, and lowest-energy-conformer `%Vbur` to
`atom_level.csv`. The Boltzmann mean uses `--temp` (298.15 K by default).

An ensemble-only run with buried volume writes `molecule_level.csv` and
`atom_level.csv`, but no bond table. When quantum outputs are also supplied,
the ensemble descriptors are merged into their normal molecule and atom
tables, and the normal bond table is retained. No raw conformer table or
separate ensemble results folder is created. Buried volumes use DBSTEP's
Bondi-radius table, molDscript's existing 1.17 scale factor, a 0.25 Angstrom
grid, and non-hydrogen occupancy atoms by default.

## Core Inputs & Flags
- `--opt PATH` - baseline optimization files and conformer metadata for quantum-output workflows; not required for single-point-only or ensemble-only analysis.
- `--spc PATH` - single-point energies that replace optimization SCF energies.
- `--nbo`, `--nmr`, `--charges`, `--fmo` PATH - add module-specific descriptors; pair with `--suffix_*` to specify filename tokens specific to calculation type (required for proper comformer matching). 
- `--fukui_neutral`, `--fukui_reduced`, `--fukui_oxidized` PATH - supply all three charge states for vertical IE/EA and condensed Fukui functions. Again, pair with `--suffix_*` for proper conformer matching.
- `--substructure SMARTS` - limit atom/bond descriptors to a SMARTS match; combine with `--volume` or `--vall` and optional `--radius` list for DBSTEP buried volumes.
- `--ensemble PATH` - analyze a multi-frame `.xyz` file, or recursively find standard CREST and GOAT ensemble files below a directory. This can be used alone or alongside quantum-output modules.
- `--suffix_ensemble TAG` - remove a trailing filename tag when matching an ensemble to existing quantum-output molecule keys.
- `--ensemble_radii "[3.5]"` - sphere radii (Angstrom) for atom-level buried-volume summaries; every atom is evaluated as the center. Set `[]` to omit this comparatively expensive descriptor.
- `--ensemble_grid 0.25` - buried-volume voxel spacing in Angstrom.
- `--ensemble_include_h` - allow hydrogen atoms to contribute to buried-volume occupancy. Hydrogen atoms still receive centered buried-volume descriptors when this flag is omitted.
- `--ensemble_exclude "[1, 2]"` - additional 1-based atoms that cannot contribute to buried-volume occupancy. Excluded atoms still receive their own centered descriptors.

Per-atom ensemble descriptors require the same atom identity and ordering in
every XYZ frame, including among atoms of the same element. Changes in element
order are rejected. A permutation of two same-element atoms is usually
undetectable from XYZ element labels alone, but remains unsupported because it
mixes atom identities across conformers.
- `--boltz`, `--min_max`, `--lowe` - compute Boltzmann-weighted averages, min/max/range tables (using `--cut`), and lowest-energy snapshots for quantum-chemistry conformer rows. These reducers are ignored for ensemble-only XYZ input because its steric ranges are already summarized internally. Adjust `--temp` (K) as needed.
- `--temp FLOAT` - temperature in Kelvin for standard conformer reducers and the atom-level ensemble `%Vbur` Boltzmann mean.
- `--output PREFIX` - prepend every generated filename; append a slash to target a directory. Use `--no_mol`, `--no_atom`, `--no_bond`, or `--no_bond_filter` to tailor CSV output.
- `--workers N` - parse independent quantum output files in parallel. Start with a modest value such as `--workers 4` for large batches, then increase if memory use is acceptable.
- `--write_args arguments.txt` - save the effective options for the current run as a reusable `--varfile`.

## Output Artefacts
- `molecule_level.csv`, `bond_level.csv`, `atom_level.csv` - aligned descriptors per calculation, bond pair, or atom.
- `ensemble_*.csv`, `boltzmann_weights.csv` - created when `--boltz` is enabled.
- `min_max_range_*.csv`, `lowest_energy_*.csv` - created when `--min_max` or `--lowe` are requested.
- `--ensemble` adds radius-of-gyration and shape min/max/range columns to `molecule_level.csv`, and per-atom buried-volume min/max/Boltzmann-mean/lowest-energy columns to `atom_level.csv` when radii are requested. It does not create additional ensemble-specific CSV files or directories.
- `MOLDSCRIPT.dat` - a single run log capturing provenance, parsed files, module sections, and CPU-time summaries.

## Documentation
The Read the Docs site (coming soon) will provide the full user guide: [https://moldscript.readthedocs.io](https://moldscript.readthedocs.io)

## Dependencies
Key Python dependencies include `pandas`, `cclib` (latest GitHub version for the most up-to-date package compatability), `dbstep`, `rdkit`, `networkx`, `numpy`, `periodictable`, `rich`, and `tqdm`.

## Supported Quantum Packages
- Gaussian
- ORCA
- xTB (optimizations)

## Testing
Run `pytest -v` from the project root to execute the test suite.

## Acknowledgements
This work was carried out in the [Paton Laboratory at Colorado State University](https://patonlab.colostate.edu), supported by the [NSF Center for Computer-Assisted Synthesis](https://ccas.nd.edu/) (grant [CHE-1925607](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2202693&HistoricalAwards=false)).

Contributors include [Shree Sowndarya](https://github.com/shreesowndarya), [Jake King](https://github.com/j77king), and [Robert Paton](https://github.com/bobbypaton).
