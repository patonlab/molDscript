"""End-to-end CLI smoke tests.

These tests invoke `python -m moldscript` as a subprocess and check that the
expected CSVs land on disk with the expected schema. The point isn't numerical
fidelity (test_moldscript.py covers that) — it's catching wiring regressions
in argument parsing, module orchestration, and CSV emission.
"""
import os
import subprocess
import sys

import pandas as pd


REPO_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), ".."))
ARBR_OPT = os.path.join("moldscript", "examples", "arbr", "opt")
ARBR_NMR = os.path.join("moldscript", "examples", "arbr", "nmr")
ARBR_OPT_ABS = os.path.join(REPO_ROOT, ARBR_OPT)
ARBR_NMR_ABS = os.path.join(REPO_ROOT, ARBR_NMR)


def run_cli(args, cwd=REPO_ROOT):
    """Run `python -m moldscript <args>` and return CompletedProcess."""
    return subprocess.run(
        [sys.executable, "-m", "moldscript", *args],
        cwd=cwd,
        capture_output=True,
        text=True,
        timeout=120,
    )


def test_cli_help_exits_zero():
    result = run_cli(["--help"])
    assert result.returncode == 0
    assert "MOLDSCRIPT" in result.stdout
    assert "--opt" in result.stdout
    assert "--suffix_opt" in result.stdout


def test_cli_opt_plus_nmr_produces_csvs(tmp_path):
    prefix = str(tmp_path) + os.sep
    result = run_cli([
        "--opt", ARBR_OPT,
        "--nmr", ARBR_NMR,
        "--suffix_nmr", "nmr",
        "--output", prefix,
    ])
    assert result.returncode == 0, f"CLI failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    for csv in ("molecule_level.csv", "bond_level.csv", "atom_level.csv"):
        path = tmp_path / csv
        assert path.exists(), f"expected {csv} to be created"
        assert path.stat().st_size > 0

    # Spot-check schemas — these are the public-API columns downstream ML
    # consumers depend on.
    mol = pd.read_csv(tmp_path / "molecule_level.csv")
    assert set(["filename", "smiles"]).issubset(mol.columns)
    assert (mol["filename"] == "arbr31_wb97xd").any()

    atom = pd.read_csv(tmp_path / "atom_level.csv")
    assert {"filename", "atom_index", "atom_type", "nmr_shielding"}.issubset(atom.columns)

    bond = pd.read_csv(tmp_path / "bond_level.csv")
    assert {"filename", "atom1_idx", "atom2_idx", "bond_length"}.issubset(bond.columns)


def test_cli_writes_module_audit_logs(tmp_path):
    env = {**os.environ, "PYTHONPATH": REPO_ROOT}
    result = subprocess.run(
        [sys.executable, "-m", "moldscript",
         "--opt", ARBR_OPT_ABS,
         "--nmr", ARBR_NMR_ABS,
         "--suffix_nmr", "nmr"],
        cwd=str(tmp_path),
        env=env,
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert result.returncode == 0, f"CLI failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    # The Logger contract: each module emits MOLDSCRIPT_<MODULE>.dat with a
    # non-empty audit trail. This test guards against future refactors that
    # accidentally bypass log.write.
    for module in ("OPT", "NMR"):
        log = tmp_path / f"MOLDSCRIPT_{module}.dat"
        assert log.exists(), f"expected {log.name} to be created"
        contents = log.read_text()
        assert "MOLDSCRIPT" in contents
        assert "Command line used in MOLDSCRIPT" in contents


def test_cli_output_prefix_reaches_dat_files(tmp_path):
    """--output should control both the CSV prefix and the per-module .dat
    audit logs. Run with cwd ≠ --output to prove the .dat tracks --output
    and not cwd."""
    cwd_dir = tmp_path / "cwd"
    out_dir = tmp_path / "out"
    cwd_dir.mkdir()
    out_dir.mkdir()
    out_prefix = str(out_dir) + os.sep
    env = {**os.environ, "PYTHONPATH": REPO_ROOT}
    result = subprocess.run(
        [sys.executable, "-m", "moldscript",
         "--opt", ARBR_OPT_ABS,
         "--output", out_prefix],
        cwd=str(cwd_dir),
        env=env,
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert result.returncode == 0, f"CLI failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    # CSVs and .dat must both land at the --output prefix, not in cwd.
    assert (out_dir / "molecule_level.csv").exists()
    assert (out_dir / "MOLDSCRIPT_OPT.dat").exists()
    assert not (cwd_dir / "MOLDSCRIPT_OPT.dat").exists(), (
        "MOLDSCRIPT_OPT.dat leaked into cwd instead of the --output prefix"
    )


def test_cli_rejects_unknown_flag():
    result = run_cli(["--definitely_not_a_flag", "foo"])
    assert result.returncode != 0
    assert "unrecognized arguments" in result.stderr or "unrecognized" in result.stderr


def test_cli_varfile_loads_defaults(tmp_path):
    varfile = tmp_path / "inputs.txt"
    varfile.write_text(
        "opt : " + ARBR_OPT + "\n"
        "nmr : " + ARBR_NMR + "\n"
        "suffix_nmr : nmr\n"
        "output : " + str(tmp_path) + os.sep + "\n"
    )
    result = run_cli(["--varfile", str(varfile)])
    assert result.returncode == 0, f"CLI failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
    assert (tmp_path / "molecule_level.csv").exists()


def test_cli_flag_overrides_varfile(tmp_path):
    """CLI --output should win over an output key in the varfile."""
    bogus_prefix = str(tmp_path / "WRONG_") + os.sep
    real_prefix = str(tmp_path) + os.sep
    varfile = tmp_path / "inputs.txt"
    varfile.write_text(
        "opt : " + ARBR_OPT + "\n"
        "output : " + bogus_prefix + "\n"
    )
    result = run_cli([
        "--varfile", str(varfile),
        "--output", real_prefix,
    ])
    assert result.returncode == 0
    assert (tmp_path / "molecule_level.csv").exists()
    assert not (tmp_path / "WRONGmolecule_level.csv").exists()
