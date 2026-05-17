"""Tests for the conformer ensemble modules: boltz, min_max, lowe.

These modules read mol/atom/bond CSVs already emitted by get_df and
produce per-molecule aggregated CSVs. They have no test coverage upstream,
even though they own the most numerical logic in the codebase. The
fixture below is a hand-crafted two-conformer system so the expected
output is analytically verifiable.
"""
import os
from pathlib import Path

import pandas as pd
import pytest

from moldscript.boltz import boltz
from moldscript.lowe import lowe
from moldscript.min_max import min_max


# Two conformers of a hypothetical "mol1". The 0.001 Ha energy gap gives
# Boltzmann weights of ~0.257 / 0.743 at 298.15 K — both above the
# default --cut threshold so min_max keeps both, but distinct enough that
# weighting is non-trivial.
E_LOW = -100.001  # Hartree (conf2, more stable)
E_HIGH = -100.0   # Hartree (conf1, less stable)
W_HIGH = 0.257479  # Boltzmann weight of conf1 (less stable)
W_LOW = 0.742521   # Boltzmann weight of conf2 (more stable)


@pytest.fixture
def ensemble_inputs(tmp_path):
    """Write a two-conformer fixture (mol1_conf1, mol1_conf2) to tmp_path
    and return (prefix, energies_df). Mirrors the schema get_df emits."""
    prefix = str(tmp_path) + os.sep

    pd.DataFrame({
        "filename": ["mol1_conf1", "mol1_conf2"],
        "smiles":   ["C(=O)O",     "C(=O)O"],
        "dipole":   [2.0,          3.0],
    }).to_csv(tmp_path / "molecule_level.csv", index=False)

    pd.DataFrame({
        "filename":   ["mol1_conf1", "mol1_conf1", "mol1_conf2", "mol1_conf2"],
        "atom_index": [1,            2,            1,            2],
        "atom_type":  ["C",          "O",          "C",          "O"],
        "charge":     [0.10,        -0.20,         0.15,        -0.30],
    }).to_csv(tmp_path / "atom_level.csv", index=False)

    pd.DataFrame({
        "filename":    ["mol1_conf1",  "mol1_conf2"],
        "atom1_idx":   [1,             1],
        "atom1":       ["C",           "C"],
        "atom2_idx":   [2,             2],
        "atom2":       ["O",           "O"],
        "bond_length": [1.21,          1.22],
    }).to_csv(tmp_path / "bond_level.csv", index=False)

    energies = pd.DataFrame({
        "filename":  ["mol1_conf1", "mol1_conf2"],
        "scfenergy": [E_HIGH,       E_LOW],
    })
    return prefix, energies, tmp_path


def test_boltz_molecule_level_weighted_average(ensemble_inputs):
    prefix, energies, tmp_path = ensemble_inputs
    boltz(temp=298.15, prefix=prefix, energies=energies)

    out = pd.read_csv(tmp_path / "ensemble_molecule_level.csv")
    assert list(out["filename"]) == ["mol1"]
    # Weighted dipole = 0.2575*2.0 + 0.7425*3.0 = 2.7425
    assert out["dipole"].iloc[0] == pytest.approx(2.7425, abs=1e-4)


def test_boltz_writes_weights_csv(ensemble_inputs):
    prefix, energies, tmp_path = ensemble_inputs
    boltz(temp=298.15, prefix=prefix, energies=energies)

    weights = pd.read_csv(tmp_path / "boltzmann_weights.csv")
    weights = weights.set_index("filename")["boltzmann_weight"]
    # Higher-energy conformer gets less weight; weights sum to 1.
    assert weights["mol1_conf1"] == pytest.approx(W_HIGH, abs=1e-4)
    assert weights["mol1_conf2"] == pytest.approx(W_LOW, abs=1e-4)
    assert weights.sum() == pytest.approx(1.0, abs=1e-6)


def test_boltz_atom_level_weighted_average(ensemble_inputs):
    prefix, energies, tmp_path = ensemble_inputs
    boltz(temp=298.15, prefix=prefix, energies=energies)

    out = pd.read_csv(tmp_path / "ensemble_atom_level.csv")
    by_idx = out.set_index("atom_index")
    # atom 1 (C): weighted charge = 0.2575*0.10 + 0.7425*0.15 = 0.1371
    assert by_idx.loc[1, "charge"] == pytest.approx(0.1371, abs=1e-4)
    # atom 2 (O): weighted charge = 0.2575*-0.20 + 0.7425*-0.30 = -0.2743
    assert by_idx.loc[2, "charge"] == pytest.approx(-0.2743, abs=1e-4)


def test_boltz_bond_level_weighted_average(ensemble_inputs):
    prefix, energies, tmp_path = ensemble_inputs
    boltz(temp=298.15, prefix=prefix, energies=energies)

    out = pd.read_csv(tmp_path / "ensemble_bond_level.csv")
    # Single bond row (C–O); weighted length = 0.2575*1.21 + 0.7425*1.22 = 1.2174
    assert len(out) == 1
    assert out["bond_length"].iloc[0] == pytest.approx(1.2174, abs=1e-4)


def test_lowe_picks_lowest_energy_conformer(ensemble_inputs):
    prefix, energies, tmp_path = ensemble_inputs
    lowe(prefix=prefix, energies=energies)

    mol = pd.read_csv(tmp_path / "lowest_energy_molecule_level.csv")
    assert list(mol["filename"]) == ["mol1"]
    # conf2 was the lower-energy structure → its dipole=3.0 wins
    assert mol["dipole"].iloc[0] == pytest.approx(3.0, abs=1e-6)

    atoms = pd.read_csv(tmp_path / "lowest_energy_atom_level.csv")
    by_idx = atoms.set_index("atom_index")
    assert by_idx.loc[1, "charge"] == pytest.approx(0.15, abs=1e-6)
    assert by_idx.loc[2, "charge"] == pytest.approx(-0.30, abs=1e-6)

    bonds = pd.read_csv(tmp_path / "lowest_energy_bond_level.csv")
    assert bonds["bond_length"].iloc[0] == pytest.approx(1.22, abs=1e-6)


def test_min_max_keeps_both_conformers_above_threshold(ensemble_inputs):
    prefix, energies, tmp_path = ensemble_inputs
    # cut=0.95 → threshold = 0.05. Both weights (0.257, 0.743) are above
    # that, so the min/max range covers both conformer values.
    min_max(temp=298.15, cut=0.95, prefix=prefix, energies=energies)

    mol = pd.read_csv(tmp_path / "min_max_range_molecule_level.csv")
    assert mol["dipole_min"].iloc[0] == pytest.approx(2.0, abs=1e-6)
    assert mol["dipole_max"].iloc[0] == pytest.approx(3.0, abs=1e-6)
    assert mol["dipole_range"].iloc[0] == pytest.approx(1.0, abs=1e-6)


def test_min_max_atom_and_bond_level(ensemble_inputs):
    prefix, energies, tmp_path = ensemble_inputs
    min_max(temp=298.15, cut=0.95, prefix=prefix, energies=energies)

    atoms = pd.read_csv(tmp_path / "min_max_range_atom_level.csv")
    by_idx = atoms.set_index("atom_index")
    assert by_idx.loc[1, "charge_min"] == pytest.approx(0.10, abs=1e-6)
    assert by_idx.loc[1, "charge_max"] == pytest.approx(0.15, abs=1e-6)
    assert by_idx.loc[1, "charge_range"] == pytest.approx(0.05, abs=1e-6)
    # Atom 2 charge ranges from -0.30 to -0.20.
    assert by_idx.loc[2, "charge_min"] == pytest.approx(-0.30, abs=1e-6)
    assert by_idx.loc[2, "charge_max"] == pytest.approx(-0.20, abs=1e-6)
    assert by_idx.loc[2, "charge_range"] == pytest.approx(0.10, abs=1e-6)

    bonds = pd.read_csv(tmp_path / "min_max_range_bond_level.csv")
    assert bonds["bond_length_min"].iloc[0] == pytest.approx(1.21, abs=1e-6)
    assert bonds["bond_length_max"].iloc[0] == pytest.approx(1.22, abs=1e-6)
    assert bonds["bond_length_range"].iloc[0] == pytest.approx(0.01, abs=1e-6)


def test_min_max_filters_conformer_below_threshold(tmp_path):
    """With a large energy gap, the high-energy conformer's Boltzmann
    weight falls below the cut threshold and should be excluded from the
    min/max range."""
    prefix = str(tmp_path) + os.sep

    # 0.01 Ha gap → weight ratio ~exp(-10.59) ~ 2.5e-5; conf1 gets filtered
    # at default cut=0.95 (threshold 0.05).
    pd.DataFrame({
        "filename": ["mol1_conf1", "mol1_conf2"],
        "smiles":   ["C(=O)O",     "C(=O)O"],
        "dipole":   [2.0,          3.0],
    }).to_csv(tmp_path / "molecule_level.csv", index=False)

    pd.DataFrame({
        "filename":   ["mol1_conf1", "mol1_conf2"],
        "atom_index": [1,            1],
        "atom_type":  ["C",          "C"],
        "charge":     [0.10,         0.15],
    }).to_csv(tmp_path / "atom_level.csv", index=False)

    pd.DataFrame({
        "filename":    ["mol1_conf1",  "mol1_conf2"],
        "atom1_idx":   [1,             1],
        "atom1":       ["C",           "C"],
        "atom2_idx":   [2,             2],
        "atom2":       ["O",           "O"],
        "bond_length": [1.21,          1.22],
    }).to_csv(tmp_path / "bond_level.csv", index=False)

    energies = pd.DataFrame({
        "filename":  ["mol1_conf1", "mol1_conf2"],
        "scfenergy": [-100.0,       -100.01],
    })
    min_max(temp=298.15, cut=0.95, prefix=prefix, energies=energies)

    mol = pd.read_csv(tmp_path / "min_max_range_molecule_level.csv")
    # Only conf2 (dipole=3.0) survives the threshold → min == max, range == 0.
    assert mol["dipole_min"].iloc[0] == pytest.approx(3.0, abs=1e-6)
    assert mol["dipole_max"].iloc[0] == pytest.approx(3.0, abs=1e-6)
    assert mol["dipole_range"].iloc[0] == pytest.approx(0.0, abs=1e-6)
