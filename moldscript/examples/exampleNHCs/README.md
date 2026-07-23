# NHC example

Run this example from the `exampleNHCs` directory:

```shell
python -m moldscript \
  --opt optimisation --suffix_opt optimisation \
  --spc singlepoint --suffix_spc singlepoint \
  --charges singlepoint --suffix_charges singlepoint \
  --fmo singlepoint --suffix_fmo singlepoint \
  --nbo singlepoint --suffix_nbo singlepoint \
  --nmr singlepoint --suffix_nmr singlepoint \
  --output results/
```

The ASE optimization logs contain optimizer energies but no coordinates, so
molDscript initializes structures from the matching Gaussian single-point
files. Those single-point files also provide the SPC energy, Mulliken charges,
frontier orbitals and moments, natural charges, and NMR shielding tensors.

The `crest_conformers` directory is not used by this example.
