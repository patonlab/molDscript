# NHC example

Run this example from the `exampleNHCs` directory:

```shell
python -m moldscript \
  --spc singlepoint --suffix_spc singlepoint \
  --charges singlepoint --suffix_charges singlepoint \
  --fmo singlepoint --suffix_fmo singlepoint \
  --nbo singlepoint --suffix_nbo singlepoint \
  --nmr singlepoint --suffix_nmr singlepoint \
  --ensemble crest_conformers \
  --ensemble_radii "[3.5]"
```

The Gaussian single-point files initialize the molecular structures and provide
the SPC energy, frontier orbitals and moments, natural charges, and NMR
shielding tensors. The multi-frame XYZ files add radius-of-gyration and shape
min/max/range columns to the same `molecule_level.csv`. For every atom, a
3.5-Angstrom sphere adds buried-volume minimum, maximum, Boltzmann mean, and
lowest-energy-conformer columns to `atom_level.csv`.

The XYZ ensembles can also be analyzed without the single-point files:

```shell
python -m moldscript \
  --ensemble crest_conformers \
  --ensemble_radii "[3.5]"
```

The combined run writes the normal molecule, atom, and bond tables, with the
ensemble values integrated into the first two. The XYZ-only form writes
`molecule_level.csv` and `atom_level.csv`, but no bond table. It creates no
separate ensemble results folder or raw conformer table. Buried volume is
atom-level only; the molecule table contains min, max, and range for radius of
gyration and shape anisotropy.

Every atom, not only Ni, is evaluated as a buried-volume center. By default,
hydrogens do not contribute to occupancy, although hydrogen atoms still
receive centered values. Use `--ensemble_include_h` to include hydrogen
occupancy, `--ensemble_exclude "[1, 2]"` to omit additional occupancy atoms,
and `--temp` to change the temperature used for the buried-volume Boltzmann
mean. All conformers must retain the same atom identity and ordering, including
same-element atoms. The default voxel spacing is 0.25 Angstrom. The
`optimisation` directory is not used for this dataset.
