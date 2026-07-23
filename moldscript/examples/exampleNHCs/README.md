# NHC example

Run this example from the `exampleNHCs` directory:

```shell
python -m moldscript \
  --spc singlepoint --suffix_spc singlepoint \
  --charges singlepoint --suffix_charges singlepoint \
  --fmo singlepoint --suffix_fmo singlepoint \
  --nbo singlepoint --suffix_nbo singlepoint \
  --nmr singlepoint --suffix_nmr singlepoint \
  --output results/
```

The Gaussian single-point files initialize the molecular structures and provide
the SPC energy, frontier orbitals and moments, natural charges, and NMR
shielding tensors.

The `optimisation` and `crest_conformers` directories are not used by this
example.
