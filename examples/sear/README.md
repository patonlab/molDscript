# SEAr regioselectivity benchmark

End-to-end example: use moldscript's Fukui descriptors (computed from
Gaussian DFT single points on neutral, +1, and -1 states) to rank the
aromatic CH positions of substituted arenes and heteroarenes, and check
whether the top-ranked carbon matches the experimentally observed major
product of electrophilic aromatic substitution.

The 18-substrate demo set is hand-curated from textbook EAS chemistry
(activated, deactivated, and halogen-substituted benzenes plus five-membered
heteroarenes). It is **not** the full Kromann/Jensen RegioSQM benchmark —
that's 525 substrates in `SI_RegioML.zip` (150 MB) on the Jensen group's
ERDA share — but the harness is dataset-agnostic: drop additional rows into
[substrates.csv](substrates.csv) and the rest of the pipeline scales.

## Workflow

```
substrates.csv                                   # SMILES + observed-site tags (atom-mapped)
        |
        v
prepare_opt.py    (RDKit / ETKDG / MMFF)         # writes gaussian/opt/*.com + targets.json
        |
        v
[run Gaussian opt jobs externally]               # produces gaussian/opt/*.log
        |
        v
prepare_singlepoints.py  (cclib reads opt geom)  # writes gaussian/{neutral,oxidized,reduced}/*.com
        |
        v
[run Gaussian SP jobs externally]                # 3x per substrate
        |
        v
python -m moldscript --varfile arguments.txt     # writes results/sear_atom_level.csv (+ mol, bond)
        |
        v
analyze.py                                       # top-1/top-2 accuracy by descriptor & class
```

## Running it

```shell
cd examples/sear

# 1. Generate opt inputs and resolve observed atom indices from the SMARTS tags.
python prepare_opt.py
#   anisole              observed_atom_index=[6]  -> gaussian/opt/anisole_conf1.com
#   ...

# 2. Submit gaussian/opt/*.com to your scheduler. Wait for *.log files.

# 3. Generate the three Fukui single points at the optimized geometry.
python prepare_singlepoints.py

# 4. Submit gaussian/neutral/*.com, gaussian/oxidized/*.com, gaussian/reduced/*.com.

# 5. Aggregate descriptors with moldscript.
python -m moldscript --varfile arguments.txt

# 6. Score predictions.
python analyze.py
```

`analyze.py` ranks the aromatic carbons of each substrate by:

- `fplus`   — electrophilic Fukui f<sup>+</sup> = −[q(N+1) − q(N)] (predicts attack by electrophiles)
- `frad`    — radical Fukui f<sup>0</sup> = (f<sup>+</sup> + f<sup>−</sup>) / 2
- `dual`    — dual descriptor f<sup>2</sup> = f<sup>+</sup> − f<sup>−</sup> (Morell)

and reports top-1 / top-2 hit rate vs. the experimentally observed site for
each descriptor, broken down by substrate class.

## Atom indexing convention

RDKit's `AddHs()` appends hydrogens **after** heavy atoms, so heavy-atom
indices survive the SMILES → 3D → Gaussian → cclib roundtrip unchanged. The
1-indexed `atom_index` column in `sear_atom_level.csv` therefore equals
`rdkit_heavy_atom_idx + 1` for any non-hydrogen atom. `targets.json` stores
the observed sites in that same 1-indexed convention.

The `observed_smarts` column in `substrates.csv` uses atom-mapped SMARTS
(`[cH:1]`) to mark the major-product carbon. `prepare_opt.py` resolves the
map number to the corresponding heavy-atom index via RDKit substructure
matching; molecules with rotational symmetry (e.g. anisole) match the
SMARTS multiple times pointing to equivalent positions, all of which count
as correct in `analyze.py`.

## Scaling to the full RegioSQM / RegioML benchmark

The 525-substrate benchmark with experimental site labels lives in
`SI_RegioML.zip` at
[sid.erda.dk/sharelink/HypB1igzDl](https://sid.erda.dk/sharelink/HypB1igzDl)
(Jensen group). After downloading and extracting:

1. Convert the dataset's site labels into the `observed_smarts` convention
   used here (one atom-mapped SMARTS per substrate, tagged at the major-product C).
2. Concatenate onto `substrates.csv`.
3. Re-run the same six-step workflow above.

Compute budget scales as ~4 single-point jobs × N substrates × N conformers.
For the demo set (N=18, 1 conformer each) at B3LYP/6-31G(d) this is ~72
small jobs (~minutes each on a workstation). For 525 substrates with
modest conformer sampling, plan for HPC.

## References

- Kromann, J. C., Jensen, J. H., Kruszyk, M., Jessing, M., & Jorgensen, M. (2018).
  Fast and accurate prediction of the regioselectivity of electrophilic aromatic
  substitution reactions. *Chem. Sci.*, 9(3), 660-665.
  https://doi.org/10.1039/C7SC04156J
- Ree, N., Goller, A. H., & Jensen, J. H. (2022). RegioML: predicting the
  regioselectivity of electrophilic aromatic substitution reactions using
  machine learning. *Digital Discovery*, 1(2), 108-114.
  https://doi.org/10.1039/D1DD00032B
- Morell, C., Grand, A., & Toro-Labbe, A. (2005). New dual descriptor for
  chemical reactivity. *J. Phys. Chem. A*, 109(1), 205-212.
  https://doi.org/10.1021/jp046577a
- Wang, B., Geerlings, P., & Heidar-Zadeh, F. (2025). Exploring intrinsic bond
  properties with the Fukui matrix from conceptual density matrix functional
  theory. *J. Chem. Theory Comput.* https://doi.org/10.1021/acs.jctc.4c01627
