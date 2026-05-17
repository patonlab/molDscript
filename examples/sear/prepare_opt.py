"""Generate Gaussian opt inputs for the SEAr regioselectivity benchmark.

Reads `substrates.csv`, embeds each SMILES with RDKit/ETKDG, force-field
preoptimizes, and writes one Gaussian opt input per substrate to
`gaussian/opt/<name>_conf1.com`. Also resolves the experimentally observed
major-product aromatic CH(s) from the atom-mapped SMARTS column and writes
1-indexed atom indices (matching moldscript's `atom_index` convention) to
`targets.json` for `analyze.py` to consume.

The RDKit heavy-atom order is preserved by AddHs (Hs are appended), so
the cclib-parsed `atom_index` in the moldscript CSVs equals
`rdkit_heavy_atom_idx + 1` for any non-hydrogen atom.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem


HERE = Path(__file__).resolve().parent
DEFAULT_CSV = HERE / "substrates.csv"
DEFAULT_OUT = HERE / "gaussian" / "opt"
DEFAULT_TARGETS = HERE / "targets.json"

OPT_ROUTE = "#P opt B3LYP/6-31G(d) integral=ultrafine"


def embed_3d(smiles: str, seed: int = 0xC057) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles!r}")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise RuntimeError(f"ETKDG embedding failed for {smiles!r}")
    AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    return mol


def write_gaussian_opt(path: Path, mol: Chem.Mol, title: str) -> None:
    conf = mol.GetConformer()
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write("%mem=4GB\n")
        f.write("%nprocshared=4\n")
        f.write(f"%chk={path.stem}.chk\n")
        f.write(f"{OPT_ROUTE}\n\n")
        f.write(f"{title}\n\n")
        f.write("0 1\n")
        for atom in mol.GetAtoms():
            p = conf.GetAtomPosition(atom.GetIdx())
            f.write(f"{atom.GetSymbol():>2}  {p.x:>12.6f}  {p.y:>12.6f}  {p.z:>12.6f}\n")
        f.write("\n")


def observed_indices(smiles: str, observed_smarts: str) -> list[int]:
    """Return 1-indexed heavy-atom indices that the atom-map-1 query tags."""
    mol = Chem.MolFromSmiles(smiles)
    query = Chem.MolFromSmarts(observed_smarts)
    if mol is None or query is None:
        raise ValueError(f"Bad SMILES/SMARTS pair: {smiles!r}, {observed_smarts!r}")
    target = next((a.GetIdx() for a in query.GetAtoms() if a.GetAtomMapNum() == 1), None)
    if target is None:
        raise ValueError(f"No [atom:1] tag in observed_smarts: {observed_smarts!r}")
    # uniquify=False so symmetry-equivalent positions (e.g. both alpha Cs of
    # pyrrole, both rings of biphenyl) are all returned as valid observed sites.
    matches = mol.GetSubstructMatches(query, uniquify=False)
    if not matches:
        raise ValueError(f"observed_smarts did not match {smiles!r}: {observed_smarts!r}")
    return sorted({m[target] + 1 for m in matches})


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--targets", type=Path, default=DEFAULT_TARGETS)
    args = parser.parse_args(argv)

    targets: dict[str, dict] = {}
    with args.csv.open() as f:
        for row in csv.DictReader(f):
            name = row["name"]
            smiles = row["smiles"]
            mol = embed_3d(smiles)
            out_com = args.out / f"{name}_conf1.com"
            write_gaussian_opt(out_com, mol, title=f"{name} opt B3LYP/6-31G(d)")
            idx = observed_indices(smiles, row["observed_smarts"])
            targets[name] = {
                "smiles": smiles,
                "observed_atom_indices": idx,
                "expected_position": row["expected_position"],
                "class": row["class"],
            }
            print(f"  {name:>20s}  observed_atom_index={idx}  -> {out_com.relative_to(HERE)}")

    args.targets.write_text(json.dumps(targets, indent=2) + "\n")
    print(f"\nWrote {len(targets)} opt inputs under {args.out.relative_to(HERE)}/")
    print(f"Wrote targets to {args.targets.relative_to(HERE)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
