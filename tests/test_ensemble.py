import datetime
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from conftest import datapath
from moldscript.argument_parser import command_line_args
from moldscript.ensemble import (
    _sphere_grid,
    boltzmann_weights,
    discover_ensemble_xyz,
    ensemble,
    parse_atom_indices,
    parse_radii,
    read_xyz_ensemble,
)
from moldscript.get_df import get_df


def _write_xyz(path, frames):
    lines = []
    for comment, atoms in frames:
        lines.extend([str(len(atoms)), str(comment)])
        lines.extend(
            f"{symbol} {x:.10f} {y:.10f} {z:.10f}"
            for symbol, x, y, z in atoms
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def test_reader_supports_crest_and_goat_comments_and_names(tmp_path):
    crest_xyz = _write_xyz(
        tmp_path / "sample_crest_conformers.xyz",
        [
            (-10.0, [("C", 0, 0, 0), ("H", 1, 0, 0)]),
            (-9.999, [("C", 0, 1, 0), ("H", 1, 1, 0)]),
        ],
    )
    goat_xyz = _write_xyz(
        tmp_path / "goat_sample.finalensemble.xyz",
        [
            (
                "-20.0 converged=true",
                [("N", 0, 0, 0), ("H", 0, 0, 1)],
            ),
            (
                "Energy = -19.999",
                [("N", 0, 1, 0), ("H", 0, 1, 1)],
            ),
        ],
    )

    crest_data = read_xyz_ensemble(crest_xyz)
    goat_data = read_xyz_ensemble(goat_xyz)

    assert crest_data.name == "sample"
    assert crest_data.symbols == ("C", "H")
    assert crest_data.coordinates.shape == (2, 2, 3)
    assert crest_data.energies_hartree.tolist() == [-10.0, -9.999]
    assert goat_data.name == "goat_sample"
    assert goat_data.energies_hartree.tolist() == [-20.0, -19.999]


def test_directory_discovery_selects_only_standard_ensemble_files(tmp_path):
    crest_file = tmp_path / "molecule_crest_conformers.xyz"
    goat_file = tmp_path / "molecule.finalensemble.xyz"
    crest_file.write_text("", encoding="utf-8")
    goat_file.write_text("", encoding="utf-8")
    (tmp_path / "crest_best.xyz").write_text("", encoding="utf-8")
    (tmp_path / "molecule.globalminimum.xyz").write_text("", encoding="utf-8")
    (tmp_path / "molecule.finalensemble.globaliter.1.xyz").write_text(
        "", encoding="utf-8"
    )

    assert discover_ensemble_xyz(tmp_path) == sorted(
        [crest_file.resolve(), goat_file.resolve()], key=str
    )


def test_cli_uses_generic_ensemble_options(monkeypatch):
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "moldscript",
            "--ensemble",
            "xyz_ensembles",
            "--ensemble_radii",
            "[3.0]",
            "--ensemble_grid",
            "0.5",
            "--suffix_ensemble",
            "GOAT",
        ],
    )

    args = command_line_args()

    assert args.ensemble == "xyz_ensembles"
    assert args.ensemble_radii == "[3.0]"
    assert args.ensemble_grid == pytest.approx(0.5)
    assert args.suffix_ensemble == "GOAT"
    assert not hasattr(args, "crest")
    assert not hasattr(args, "ensemble_center")


@pytest.mark.parametrize(
    "contents, message",
    [
        (
            "2\n-10.0\nC 0 0 0\n",
            "truncated after 1 of 2 atom rows",
        ),
        (
            "2\n-10.0\nC 0 0 0\nH 1 0 0\n"
            "2\n-9.9\nH 0 0 0\nC 1 0 0\n",
            "changes atom ordering",
        ),
        (
            "1\nnot-an-energy\nC 0 0 0\n",
            "expected a Hartree energy",
        ),
        (
            "1\n-10.0\nC 0 0 0\n"
            "2\n-9.9\nC 0 0 0\nH 1 0 0\n",
            "frame 2 contains 2 atoms",
        ),
        (
            "1\n-10.0\nC 0 nan 0\n",
            "coordinates must be finite",
        ),
        (
            "1\n-10.0 converged=false\nC 0 0 0\n",
            "marked converged=false",
        ),
        (
            "1\n-10.0\nXx 0 0 0\n",
            "unknown element symbol",
        ),
    ],
)
def test_reader_rejects_malformed_ensembles(tmp_path, contents, message):
    xyz = tmp_path / "bad.xyz"
    xyz.write_text(contents, encoding="utf-8")

    with pytest.raises(ValueError, match=message):
        read_xyz_ensemble(xyz)


def test_ensemble_numeric_options_are_validated():
    assert parse_radii("[3.5, 5]") == [3.5, 5.0]
    assert parse_radii("[]") == []
    assert parse_atom_indices("[1, 3]") == [0, 2]
    with pytest.raises(ValueError, match="whole numbers"):
        parse_atom_indices([1.5])
    with pytest.raises(ValueError, match="positive"):
        parse_radii([0])


def test_boltzmann_weights_are_normalized_shift_invariant_and_temperature_aware():
    energies = np.asarray([-10.0, -9.999, -9.998])
    weights = boltzmann_weights(energies, temperature=298.15)
    shifted = boltzmann_weights(energies + 123.4, temperature=298.15)
    hotter = boltzmann_weights(energies, temperature=1000.0)

    assert weights.sum() == pytest.approx(1.0)
    assert np.allclose(weights, shifted)
    assert weights[0] > weights[1] > weights[2]
    assert hotter[0] < weights[0]
    with pytest.raises(ValueError, match="positive"):
        boltzmann_weights(energies, temperature=0)


def test_buried_volume_grid_is_centered_and_radius_independent():
    small_points, _ = _sphere_grid(1.0, 0.3)
    large_points, large_squared_radius = _sphere_grid(1.37, 0.3)
    selected_large_points = large_points[
        large_squared_radius <= 1.0 + 1.0e-12
    ]

    assert np.any(np.all(small_points == 0.0, axis=1))
    rounded_points = {
        tuple(np.round(point, 12)) for point in small_points
    }
    assert all(
        tuple(np.round(-point, 12)) in rounded_points
        for point in small_points
    )
    assert {
        tuple(np.round(point, 12)) for point in selected_large_points
    } == rounded_points
    with pytest.raises(ValueError, match="safety limit"):
        _sphere_grid(5.0, 0.01)


def test_steric_ranges_are_invariant_to_rigid_motion(tmp_path):
    first = [
        ("Ni", 0.0, 0.0, 0.0),
        ("C", 2.0, 0.0, 0.0),
        ("N", 0.0, 2.0, 0.0),
        ("H", 0.0, 0.0, 1.0),
    ]
    axis = np.asarray([1.0, 2.0, 3.0])
    axis /= np.linalg.norm(axis)
    angle = 0.731
    cross_product = np.array(
        [
            [0.0, -axis[2], axis[1]],
            [axis[2], 0.0, -axis[0]],
            [-axis[1], axis[0], 0.0],
        ]
    )
    rotation = (
        np.eye(3) * np.cos(angle)
        + (1.0 - np.cos(angle)) * np.outer(axis, axis)
        + np.sin(angle) * cross_product
    )
    transformed = (
        np.asarray([atom[1:] for atom in first]) @ rotation.T
        + np.asarray([5.0, -2.0, 1.0])
    )
    second = [
        (atom[0], *coordinates)
        for atom, coordinates in zip(first, transformed)
    ]
    xyz = _write_xyz(
        tmp_path / "rigid.xyz",
        [(-20.0, first), (-19.999, second)],
    )

    result = ensemble(
        xyz,
        radii=[3.0],
        grid=0.4,
        create_dat=False,
    )
    entry = result.file_data["rigid"]
    mol = entry["mol"]

    expected_molecule_sterics = {
        "ensemble_radius_of_gyration_min_angstrom",
        "ensemble_radius_of_gyration_max_angstrom",
        "ensemble_radius_of_gyration_range_angstrom",
        "ensemble_shape_anisotropy_min",
        "ensemble_shape_anisotropy_max",
        "ensemble_shape_anisotropy_range",
    }
    expected_atom_sterics = {
        "ensemble_buried_volume_r_3A_min_percent",
        "ensemble_buried_volume_r_3A_max_percent",
        "ensemble_buried_volume_r_3A_boltzmann_mean_percent",
        "ensemble_buried_volume_r_3A_lowest_energy_percent",
    }
    assert set(mol) == {
        "smiles",
        "scfenergy",
        *expected_molecule_sterics,
    }
    assert set(entry["atom"]) == {"atomnos", *expected_atom_sterics}
    assert mol["ensemble_radius_of_gyration_range_angstrom"] < 1.0e-9
    assert mol["ensemble_shape_anisotropy_range"] < 1.0e-9
    atom = entry["atom"]
    minimum = atom["ensemble_buried_volume_r_3A_min_percent"]
    maximum = atom["ensemble_buried_volume_r_3A_max_percent"]
    mean = atom[
        "ensemble_buried_volume_r_3A_boltzmann_mean_percent"
    ]
    lowest = atom[
        "ensemble_buried_volume_r_3A_lowest_energy_percent"
    ]
    np.testing.assert_allclose(maximum, minimum, atol=1.0e-12)
    np.testing.assert_allclose(mean, minimum, atol=1.0e-12)
    np.testing.assert_allclose(lowest, minimum, atol=1.0e-12)
    assert entry["bond"] == {}


def test_steric_ranges_are_invariant_to_global_ensemble_rotation(tmp_path):
    frames = [
        [
            ("Ni", 0.0, 0.0, 0.0),
            ("C", 2.0, 0.0, 0.0),
            ("N", 0.0, 2.0, 0.0),
            ("H", 0.0, 0.0, 1.0),
        ],
        [
            ("Ni", 0.0, 0.0, 0.0),
            ("C", 2.3, 0.1, 0.0),
            ("N", 0.2, 1.7, 0.1),
            ("H", 0.0, 0.2, 1.2),
        ],
    ]
    axis = np.asarray([2.0, -1.0, 3.0])
    axis /= np.linalg.norm(axis)
    angle = 0.613
    cross_product = np.array(
        [
            [0.0, -axis[2], axis[1]],
            [axis[2], 0.0, -axis[0]],
            [-axis[1], axis[0], 0.0],
        ]
    )
    rotation = (
        np.eye(3) * np.cos(angle)
        + (1.0 - np.cos(angle)) * np.outer(axis, axis)
        + np.sin(angle) * cross_product
    )
    translation = np.asarray([-4.0, 2.5, 1.2])
    rotated_frames = []
    for frame in frames:
        transformed = (
            np.asarray([atom[1:] for atom in frame]) @ rotation.T
            + translation
        )
        rotated_frames.append(
            [
                (atom[0], *coordinates)
                for atom, coordinates in zip(frame, transformed)
            ]
        )

    original_xyz = _write_xyz(
        tmp_path / "original.xyz",
        [(-20.0, frames[0]), (-19.999, frames[1])],
    )
    rotated_xyz = _write_xyz(
        tmp_path / "rotated.xyz",
        [(-20.0, rotated_frames[0]), (-19.999, rotated_frames[1])],
    )
    original = ensemble(
        original_xyz, radii=[3.0], grid=0.3, create_dat=False
    ).file_data["original"]
    rotated = ensemble(
        rotated_xyz, radii=[3.0], grid=0.3, create_dat=False
    ).file_data["rotated"]

    for descriptor in original["mol"]:
        if descriptor.startswith("ensemble_"):
            assert rotated["mol"][descriptor] == pytest.approx(
                original["mol"][descriptor], abs=1.0e-8
            )
    for descriptor in original["atom"]:
        if descriptor.startswith("ensemble_"):
            np.testing.assert_allclose(
                rotated["atom"][descriptor],
                original["atom"][descriptor],
                atol=1.0e-8,
            )


def test_per_atom_buried_volume_summaries_use_energy_and_temperature(
    tmp_path,
):
    energies = [-10.0, -10.002, -9.9995]
    frames = [
        [
            ("C", 0.0, 0.0, 0.0),
            ("N", 1.7, 0.2, 0.0),
            ("O", -0.2, 2.0, 0.3),
        ],
        [
            ("C", 0.0, 0.0, 0.0),
            ("N", 2.4, 0.1, 0.0),
            ("O", 0.0, 1.5, 0.5),
        ],
        [
            ("C", 0.0, 0.0, 0.0),
            ("N", 1.3, -0.4, 0.2),
            ("O", -0.3, 2.5, -0.2),
        ],
    ]
    xyz = _write_xyz(
        tmp_path / "summary.xyz",
        list(zip(energies, frames)),
    )
    result = ensemble(
        xyz,
        radii=[2.5],
        grid=0.4,
        temp=350.0,
        create_dat=False,
    ).file_data["summary"]

    conformer_values = []
    for index, (energy, frame) in enumerate(zip(energies, frames)):
        single_xyz = _write_xyz(
            tmp_path / f"single_{index}.xyz",
            [(energy, frame)],
        )
        single_atom = ensemble(
            single_xyz,
            radii=[2.5],
            grid=0.4,
            temp=350.0,
            create_dat=False,
        ).file_data[f"single_{index}"]["atom"]
        conformer_values.append(
            single_atom[
                "ensemble_buried_volume_r_2_5A_min_percent"
            ]
        )
    conformer_values = np.asarray(conformer_values)
    weights = boltzmann_weights(energies, temperature=350.0)
    atom = result["atom"]

    assert not any(
        "buried_volume" in descriptor for descriptor in result["mol"]
    )
    np.testing.assert_allclose(
        atom["ensemble_buried_volume_r_2_5A_min_percent"],
        np.min(conformer_values, axis=0),
    )
    np.testing.assert_allclose(
        atom["ensemble_buried_volume_r_2_5A_max_percent"],
        np.max(conformer_values, axis=0),
    )
    np.testing.assert_allclose(
        atom[
            "ensemble_buried_volume_r_2_5A_boltzmann_mean_percent"
        ],
        np.sum(conformer_values * weights[:, None], axis=0),
    )
    np.testing.assert_allclose(
        atom[
            "ensemble_buried_volume_r_2_5A_lowest_energy_percent"
        ],
        conformer_values[1],
    )
    assert np.any(
        np.ptp(conformer_values, axis=0) > 0
    )


def test_ensemble_only_cli_writes_molecule_and_atom_tables(tmp_path):
    xyz = _write_xyz(
        tmp_path / "standalone.xyz",
        [
            (-10.0, [("C", 0, 0, 0), ("H", 1, 0, 0)]),
            (-9.999, [("C", 0, 0, 0), ("H", 1.1, 0, 0)]),
        ],
    )
    output = tmp_path / "output"
    subprocess.run(
        [
            sys.executable,
            "-m",
            "moldscript",
            "--ensemble",
            str(xyz),
            "--ensemble_radii",
            "[3.0]",
            "--ensemble_grid",
            "0.5",
            "--output",
            f"{output}/",
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    assert (output / "molecule_level.csv").exists()
    assert (output / "atom_level.csv").exists()
    assert not (output / "bond_level.csv").exists()
    molecule_df = pd.read_csv(output / "molecule_level.csv")
    atom_df = pd.read_csv(output / "atom_level.csv")
    assert not any(
        "buried_volume" in column for column in molecule_df.columns
    )
    assert {
        column
        for column in atom_df.columns
        if "buried_volume" in column
    } == {
        "ensemble_buried_volume_r_3A_min_percent",
        "ensemble_buried_volume_r_3A_max_percent",
        "ensemble_buried_volume_r_3A_boltzmann_mean_percent",
        "ensemble_buried_volume_r_3A_lowest_energy_percent",
    }
    assert not list(output.glob("ensemble_*.csv"))
    assert not list(output.glob("crest_*.csv"))


def test_programmatic_ensemble_dataframe_writes_atom_buried_volume(tmp_path):
    xyz = _write_xyz(
        tmp_path / "api.xyz",
        [
            (-10.0, [("C", 0, 0, 0), ("H", 1, 0, 0)]),
            (-9.999, [("C", 0, 0, 0), ("H", 1.1, 0, 0)]),
        ],
    )
    result = ensemble(
        xyz, radii=[3.0], grid=0.5, create_dat=False
    )
    output = tmp_path / "api_output"
    output.mkdir()

    get_df(result.file_data, prefix=f"{output}/")

    assert sorted(path.name for path in output.iterdir()) == [
        "atom_level.csv",
        "molecule_level.csv"
    ]
    molecule_df = pd.read_csv(output / "molecule_level.csv")
    atom_df = pd.read_csv(output / "atom_level.csv")
    assert not any(
        "buried_volume" in column for column in molecule_df.columns
    )
    assert len(
        [
            column
            for column in atom_df.columns
            if "buried_volume" in column
        ]
    ) == 4


@pytest.mark.parametrize("reducer_flag", ["--boltz", "--min_max", "--lowe"])
def test_ensemble_only_cli_ignores_legacy_reducers(
    tmp_path, reducer_flag
):
    xyz = _write_xyz(
        tmp_path / "standalone.xyz",
        [
            (-10.0, [("C", 0, 0, 0), ("H", 1, 0, 0)]),
            (-9.999, [("C", 0, 0, 0), ("H", 1.1, 0, 0)]),
        ],
    )
    output = tmp_path / reducer_flag[2:]
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "moldscript",
            "--ensemble",
            str(xyz),
            "--ensemble_radii",
            "[]",
            reducer_flag,
            "--output",
            f"{output}/",
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    assert "ignored for ensemble-only XYZ input" in completed.stdout
    assert sorted(path.name for path in output.iterdir()) == [
        "MOLDSCRIPT.dat",
        "molecule_level.csv",
    ]


def test_ensemble_merges_only_sterics_into_existing_tables(tmp_path):
    xyz = _write_xyz(
        tmp_path / "sample_GOAT.xyz",
        [
            (-10.0, [("C", 0, 0, 0), ("H", 1, 0, 0)]),
            (-9.999, [("C", 0, 0, 0), ("H", 1.1, 0, 0)]),
        ],
    )
    original_bond_length = np.asarray([[0.0, 1.01], [1.01, 0.0]])
    data_dict = {
        "CPU_time": [],
        "sample": {
            "mol": {"smiles": "[H]C", "scfenergy": -11.0},
            "atom": {"atomnos": np.asarray([6, 1])},
            "bond": {"bond_length": original_bond_length.copy()},
            "CPU_time": datetime.timedelta(0),
        },
    }

    result = ensemble(
        xyz,
        data_dict=data_dict,
        suffix="GOAT",
        radii=[3.0],
        grid=0.5,
        create_dat=False,
    )
    entry = result.file_data["sample"]

    assert entry["mol"]["smiles"] == "[H]C"
    assert entry["mol"]["scfenergy"] == -11.0
    assert np.array_equal(entry["bond"]["bond_length"], original_bond_length)
    assert {
        key for key in entry["atom"] if key.startswith("ensemble_")
    } == {
        "ensemble_buried_volume_r_3A_min_percent",
        "ensemble_buried_volume_r_3A_max_percent",
        "ensemble_buried_volume_r_3A_boltzmann_mean_percent",
        "ensemble_buried_volume_r_3A_lowest_energy_percent",
    }
    assert set(entry["bond"]) == {"bond_length"}
    assert {
        key for key in entry["mol"] if key.startswith("ensemble_")
    } == {
        "ensemble_radius_of_gyration_min_angstrom",
        "ensemble_radius_of_gyration_max_angstrom",
        "ensemble_radius_of_gyration_range_angstrom",
        "ensemble_shape_anisotropy_min",
        "ensemble_shape_anisotropy_max",
        "ensemble_shape_anisotropy_range",
    }

    get_df(result.file_data, prefix=f"{tmp_path}/")
    molecule_df = pd.read_csv(tmp_path / "molecule_level.csv")
    atom_df = pd.read_csv(tmp_path / "atom_level.csv")
    bond_df = pd.read_csv(tmp_path / "bond_level.csv")

    assert "ensemble_radius_of_gyration_range_angstrom" in molecule_df
    assert not any(
        "buried_volume" in column for column in molecule_df
    )
    assert len(
        [
            column
            for column in atom_df
            if column.startswith("ensemble_buried_volume_")
        ]
    ) == 4
    assert not any(column.startswith("ensemble_") for column in bond_df)
    assert not list(tmp_path.glob("crest_*.csv"))
    assert not list(tmp_path.glob("ensemble_*.csv"))


@pytest.mark.skipif(
    not Path(datapath("exampleNHCs/crest_conformers")).exists(),
    reason="exampleNHCs ensemble files are not included in this checkout",
)
def test_example_nhcs_ensemble_counts_and_compact_sterics():
    ensemble_path = Path(datapath("exampleNHCs/crest_conformers"))
    files = discover_ensemble_xyz(ensemble_path)
    expected_counts = {
        "core_1_r_10_r_10_X_X": 6,
        "core_1_r_11_r_11_X_X": 1,
        "core_1_r_12_r_1_X_X": 198,
        "core_1_r_13_r_13_X_X": 43,
        "core_1_r_1_r_1_X_X": 25,
    }

    observed_counts = {
        parsed.name: parsed.n_conformers
        for parsed in (read_xyz_ensemble(path) for path in files)
    }
    assert observed_counts == expected_counts
    assert sum(observed_counts.values()) == 273

    result = ensemble(
        ensemble_path,
        radii=[],
        workers=1,
        create_dat=False,
    )
    for name in expected_counts:
        entry = result.file_data[name]
        assert len(
            [
                key
                for key in entry["mol"]
                if key.startswith("ensemble_")
            ]
        ) == 6
        assert set(entry["atom"]) == {"atomnos"}
        assert entry["bond"] == {}

    six_conformer_file = (
        ensemble_path
        / "core_1_r_10_r_10_X_X_crest_conformers.xyz"
    )
    steric_result = ensemble(
        six_conformer_file,
        radii=[3.5],
        workers=1,
        create_dat=False,
    )
    entry = steric_result.file_data["core_1_r_10_r_10_X_X"]
    atom = entry["atom"]
    assert not any("buried_volume" in key for key in entry["mol"])
    minimum = atom[
        "ensemble_buried_volume_r_3_5A_min_percent"
    ]
    maximum = atom[
        "ensemble_buried_volume_r_3_5A_max_percent"
    ]
    mean = atom[
        "ensemble_buried_volume_r_3_5A_boltzmann_mean_percent"
    ]
    lowest = atom[
        "ensemble_buried_volume_r_3_5A_lowest_energy_percent"
    ]
    assert minimum.shape == (46,)
    assert np.all(minimum <= mean)
    assert np.all(mean <= maximum)
    assert np.all(minimum <= lowest)
    assert np.all(lowest <= maximum)
