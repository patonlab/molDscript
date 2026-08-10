.. image:: ../../moldscript/molDscript.png
    :alt: molDscript

.. image:: https://dl.circleci.com/status-badge/img/circleci/JDvVi58JeRw4LYzfeJGsjn/T7u88FmqfkcE7c6vKveH1L/tree/main.svg?style=shield&circle-token=CCIPRJ_3nGjXb4n3dHaAo6mQ67TBk_5ce95f5de89641ed836cbe55488e9b11f28c43d3
    :target: https://dl.circleci.com/status-badge/redirect/circleci/JDvVi58JeRw4LYzfeJGsjn/T7u88FmqfkcE7c6vKveH1L/tree/main

.. image:: https://img.shields.io/badge/License-MIT-yellow.svg
    :target: https://opensource.org/licenses/MIT

molDscript User Guide
=====================

Introduction
------------
MolDscript is a Python workflow that transforms Density Functional Theory (DFT) and related quantum chemistry outputs into descriptor tables suitable for machine learning, benchmarking, and mechanistic analysis. It leverages `cclib`, `RDKit`, `DBSTEP`, and `pandas` to parse Gaussian, ORCA, and xTB (optimization) log files, align matching calculations, and write molecule-, bond-, and atom-level CSV files with consistent identifiers.

Key Capabilities
----------------
- Automate ingestion of optimization, single-point, NBO, NMR, Fukui, charge, and frontier-orbital calculations without hand editing.
- Merge descriptors across conformers and calculation types into aligned CSV datasets.
- Restrict analysis to user-defined SMARTS substructures and optionally compute DBSTEP buried volumes.
- Summarize steric variation in CREST or ORCA GOAT multi-frame XYZ ensembles without creating raw conformer tables.
- Generate ensemble statistics such as Boltzmann-weighted averages, population windows, and lowest-energy snapshots.
- Produce a single audit log (``MOLDSCRIPT.dat``) alongside descriptor files for reproducibility.

Installation
------------
1. Clone the repository: ``git clone https://github.com/patonlab/molDscript.git``.
2. (Optional) create and activate a dedicated environment.
3. Install the package from the repository root: ``pip install -e .`` (or ``pip install .`` for a standard install).

Required Python dependencies are declared in ``setup.py`` and include ``pandas>=2.0.2``, ``cclib`` (latest from GitHub), ``dbstep``, ``rdkit``, ``networkx``, ``numpy``, ``periodictable``, ``rich``, and ``tqdm``. Install RDKit and Open Babel via conda-forge when pip wheels are not available:

.. code-block:: shell

    conda install -c conda-forge rdkit openbabel

Open Babel is not strictly required but provides a fallback encoder for substructure matching when RDKit cannot build the molecule directly from a log file.

Preparing Inputs
----------------
Supported calculation software:

- Gaussian (``.log``) for all modules.
- ORCA (``.out``) for NMR, frontier orbital, charge, and Fukui analyses.
- xTB optimizations for geometry-centric data and CPU time extraction.

Organize each calculation type in its own directory. A typical layout is:

.. code-block:: text

    calculations/
      opt/
      spc/
      nbo/
      nmr/
      charges/
      fmo/
      fukui/
        neutral/
        reduced/
        oxidized/

molDscript aligns files across modules by their common stem (everything before suffixes such as ``_opt.log`` or ``_nmr.log``). If your files include extra tokens, pass the appropriate ``--suffix_*`` option (for example, ``--suffix_nmr _nmr``) so the stems match. Conformer ensembles are detected when filenames share a stem and include ``_confX`` tokens; these are required for Boltzmann, min/max/range, and lowest-energy analyses.

For Fukui indices you must supply all three charge states via ``--fukui_neutral``, ``--fukui_reduced``, and ``--fukui_oxidized``. Substructure matching and buried-volume calculations require optimized geometries (``--opt``) and a SMARTS pattern (``--substructure``).

Running molDscript
------------------
Execute molDscript from the project root (or any directory where the desired output prefix is writable):

.. code-block:: shell

    python -m moldscript --opt calculations/opt --spc calculations/spc         --nbo calculations/nbo --nmr calculations/nmr         --charges calculations/charges --fmo calculations/fmo         --fukui_neutral calculations/fukui/neutral         --fukui_reduced calculations/fukui/reduced         --fukui_oxidized calculations/fukui/oxidized         --substructure "[#6]-[#8]" --volume --radius "[3.0, 3.5]"         --boltz --min_max --lowe --output results/

All descriptor CSV files are written to the working directory; ``--output`` lets you prepend a directory or filename prefix (for example, ``--output results/`` or ``--output studyA_``).

Argument Files
--------------
Any options accepted on the command line can be stored in a simple key-value text file and supplied with ``--varfile``:

.. code-block:: text

    opt: calculations/opt
    spc: calculations/spc
    nbo: calculations/nbo
    substructure: [#6]-[#8]
    volume: true
    radius: [3.0, 3.5]
    boltz: true
    output: results/

Invoke with ``python -m moldscript --varfile inputs.txt``. Explicit command-line options override values loaded from the file.

CLI Options
-----------
The following options can be combined as needed. Paths can be absolute or relative.

Core inputs
^^^^^^^^^^^
``--opt PATH``
  Directory containing optimization log/out files. Provides the structural baseline, SCF energies, and conformer metadata when used.

``--spc PATH``
  Directory of single-point energy calculations. Overrides the SCF energies captured during the optimization step.

``--suffix_opt TEXT`` / ``--suffix_spc TEXT``
  Trailing text to strip from filenames before matching stems (for example, ``_opt`` or ``_spc``).

``--ensemble PATH``
  A multi-frame XYZ file, or a directory containing standard CREST
  ``crest_conformers.xyz`` / ``*_crest_conformers.xyz`` files or ORCA GOAT
  ``*.finalensemble.xyz`` files. Nonstandard XYZ names can be passed directly.

Property modules
^^^^^^^^^^^^^^^^
``--nbo PATH``
  Add Wiberg bond indices, bond-order matrices, and natural population charges from NBO-enabled calculations.

``--nmr PATH``
  Parse isotropic NMR shielding tensors (Gaussian or ORCA output).

``--charges PATH``
  Capture all partial charge schemes reported by the calculation (Mulliken, Hirshfeld, NPA, etc.).

``--fmo PATH``
  Extract dipole moments, HOMO/LUMO energies, HOMO-LUMO gaps, chemical potential, and global electrophilicity/nucleophilicity parameters.

``--fukui_neutral PATH`` ``--fukui_reduced PATH`` ``--fukui_oxidized PATH``
  Directories for the neutral, anionic, and cationic calculations used to compute vertical IE/EA and condensed Fukui functions. Use ``--suffix_fukui_neutral``, ``--suffix_fukui_reduced``, and ``--suffix_fukui_oxidized`` when filenames contain extra tags.

``--suffix_nbo TEXT`` ``--suffix_nmr TEXT`` ``--suffix_charges TEXT`` ``--suffix_fmo TEXT``
  Optional suffix trimming for the respective modules.

Substructure and sterics
^^^^^^^^^^^^^^^^^^^^^^^^
``--substructure SMARTS``
  Restrict atom and bond descriptors to atoms matching the SMARTS pattern. The pattern is evaluated with RDKit; Open Babel (if installed) provides a fallback encoder.

``--volume``
  Compute DBSTEP buried volumes for the matched substructure atoms. Requires ``--substructure``.

``--vall``
  Compute DBSTEP buried volumes for every atom in each molecule (ignores ``--substructure``).

``--radius VALUE`` or ``--radius "[3.0, 3.5]"``
  Probe radii (in Angstrom) passed to DBSTEP. Accepts a single float or a Python-like list string.

XYZ ensemble sterics
^^^^^^^^^^^^^^^^^^^^
``--ensemble_radii "[3.5]"``
  Sphere radii used for atom-level percent-buried-volume descriptors. Every
  atom is evaluated as a center. For each atom and radius, molDscript reports
  the minimum, maximum, Boltzmann-weighted mean, and value from the
  lowest-energy conformer. Use ``[]`` to omit buried volume while retaining
  molecule-level size and shape summaries.

``--ensemble_grid 0.25``
  Buried-volume voxel spacing in Angstrom. The default is 0.25 Angstrom.

``--ensemble_include_h`` / ``--ensemble_exclude "[1, 2]"``
  Control which atoms contribute to buried-volume occupancy.
  ``--ensemble_include_h`` includes hydrogens as occupancy atoms; otherwise
  they are omitted from occupancy. ``--ensemble_exclude`` removes additional
  1-based atoms from occupancy. These controls do not remove atom rows or
  centers: every atom, including a hydrogen or excluded atom, still receives
  centered buried-volume descriptors.

``--suffix_ensemble TEXT``
  Strip a trailing XYZ filename tag when matching ensemble data to existing
  quantum-output molecule keys.

Ensemble XYZ input adds molecule-level min/max/range descriptors only for
mass-weighted radius of gyration and relative shape anisotropy; buried volume
is not a molecule-level descriptor. Requested buried volumes are written to
``atom_level.csv`` as the min, max, Boltzmann mean, and
lowest-energy-conformer value for every atom and radius. The Boltzmann mean
uses the temperature supplied by ``--temp``.

A standalone ensemble run with buried volume writes ``molecule_level.csv`` and
``atom_level.csv`` but no bond table. With ``--ensemble_radii "[]"``, it
writes only the molecule table. When combined with quantum outputs, ensemble
values are merged into the existing molecule and atom tables and the normal
bond table remains available. No raw conformer table, separate results folder,
or other ensemble-specific CSV is created.

These atom-level summaries require stable atom identity and ordering in every
frame, including among same-element atoms. Different element ordering is
rejected. Same-element permutations are usually impossible to detect from XYZ
element labels alone, but they are unsupported because they mix atom
identities across conformers.

Ensemble statistics
^^^^^^^^^^^^^^^^^^^
These reducers operate on quantum-chemistry conformer rows and are ignored for
ensemble-only XYZ input.

``--boltz``
  Produce Boltzmann-weighted averages at the temperature specified by ``--temp`` (default 298.15 K), along with the conformer weights.

``--min_max``
  Compute min/max/range descriptors for conformers whose Boltzmann weight exceeds ``1 - --cut`` (default retention threshold 0.95).

``--lowe``
  Keep only the lowest-energy conformer for each molecule.

``--temp FLOAT``
  Temperature (in Kelvin) used for Boltzmann and min/max population analyses,
  as well as the atom-level ensemble buried-volume Boltzmann mean.

``--cut FLOAT``
  Cumulative Boltzmann weight cutoff for ``--min_max`` (expressed as the retained population fraction).

Output control
^^^^^^^^^^^^^^
``--output PREFIX``
  Prefix applied to every generated file name. Use a trailing slash to direct output to a directory.

``--no_mol`` ``--no_atom`` ``--no_bond``
  Skip writing molecule-, atom-, or bond-level descriptor tables respectively.

``--no_bond_filter``
  Disable the default filtering (``bond_order >= 0.1`` when available, otherwise bond length <= 3 Angstrom).

``--varfile PATH``
  Load arguments from a key:value text file (see the example above).

Output Files
------------
Running molDscript creates the following artefacts in the working directory (or under ``--output``):

- ``molecule_level.csv`` - one row per calculation with molecular descriptors (SCF energy, SMILES, dipole, frontier orbital metrics, etc.).
- ``bond_level.csv`` - pairwise descriptors filtered by the default bond-order/length thresholds.
- ``atom_level.csv`` - atomic descriptors including charges, Fukui indices, NMR shielding, and buried volumes when requested.
- ``boltzmann_weights.csv`` plus ``ensemble_*.csv`` tables when ``--boltz`` is enabled.
- ``min_max_range_*.csv`` tables when ``--min_max`` is enabled and ``lowest_energy_*.csv`` tables when ``--lowe`` is requested.
- ``--ensemble`` writes radius-of-gyration and shape summaries directly into
  ``molecule_level.csv`` and requested per-atom buried-volume summaries into
  ``atom_level.csv``. It creates no additional ensemble-specific table or
  directory.
- ``MOLDSCRIPT.dat`` - a single run log that documents command provenance, parsed files, module sections, and CPU-time summaries.

Each run also reports the cumulative CPU time associated with the parsed quantum chemistry jobs.

Testing
-------
Execute ``pytest -v`` from the project root to run the automated test suite.

Acknowledgements
----------------
This work was carried out in the `Paton Laboratory at Colorado State University <https://patonlab.colostate.edu>`_, supported by the `NSF Center for Computer-Assisted Synthesis <https://ccas.nd.edu>`_ (grant `CHE-1925607 <https://www.nsf.gov/awardsearch/showAward?AWD_ID=2202693&HistoricalAwards=false>`_).

Contributors include `Shree Sowndarya <https://github.com/shreesowndarya>`_, `Jake King <https://github.com/j77king>`_, and `Robert Paton <https://github.com/bobbypaton>`_.
