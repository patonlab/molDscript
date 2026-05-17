"""Score SEAr regioselectivity predictions from moldscript Fukui descriptors.

For each substrate:
  1. Read aromatic carbon atom indices from the input SMILES (via RDKit).
  2. Slice `results/sear_atom_level.csv` to those atoms only.
  3. Rank aromatic Cs by candidate descriptors (fplus, frad, dual = fplus-fminus).
  4. Compare top-1 / top-2 ranking to the experimentally observed indices
     stored in `targets.json`.

Reports per-substrate hit/miss for each descriptor and an overall accuracy
table. No model fitting — this is a sanity check that the descriptor
direction matches textbook EAS regiochemistry on the demo set.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem


HERE = Path(__file__).resolve().parent
DEFAULT_TARGETS = HERE / "targets.json"
DEFAULT_ATOM_CSV = HERE / "results" / "sear_atom_level.csv"

DESCRIPTORS = ("fplus", "frad", "dual")


def aromatic_carbon_indices(smiles: str) -> list[int]:
    """1-indexed heavy-atom indices of aromatic carbons in the SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    return [a.GetIdx() + 1 for a in mol.GetAtoms() if a.GetSymbol() == "C" and a.GetIsAromatic()]


def filename_to_name(fname: str) -> str:
    stem = Path(fname).stem
    return stem.removesuffix("_conf1")


def rank_sites(atom_df: pd.DataFrame, descriptor: str) -> list[int]:
    """Return aromatic-C atom_indices ordered by descending descriptor."""
    return atom_df.sort_values(descriptor, ascending=False)["atom_index"].astype(int).tolist()


def score(targets: dict, atom_csv: Path) -> pd.DataFrame:
    df = pd.read_csv(atom_csv)
    if "dual" not in df.columns and {"fplus", "fminus"}.issubset(df.columns):
        df["dual"] = df["fplus"] - df["fminus"]
    df["name"] = df["filename"].map(filename_to_name)

    rows = []
    for name, meta in targets.items():
        sub = df[df["name"] == name]
        if sub.empty:
            print(f"  ! {name}: no rows in atom_level.csv", file=sys.stderr)
            continue
        arom_idx = set(aromatic_carbon_indices(meta["smiles"]))
        sub = sub[sub["atom_index"].astype(int).isin(arom_idx)]
        if sub.empty:
            print(f"  ! {name}: no aromatic-C rows after filtering", file=sys.stderr)
            continue
        observed = set(meta["observed_atom_indices"])
        row = {"name": name, "expected_position": meta["expected_position"], "class": meta["class"]}
        for desc in DESCRIPTORS:
            if desc not in sub.columns:
                row[f"top1_{desc}"] = None
                row[f"top2_{desc}"] = None
                continue
            ranked = rank_sites(sub, desc)
            row[f"top1_{desc}"] = int(ranked[0] in observed) if ranked else None
            row[f"top2_{desc}"] = int(any(i in observed for i in ranked[:2])) if ranked else None
        rows.append(row)
    return pd.DataFrame(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--targets", type=Path, default=DEFAULT_TARGETS)
    parser.add_argument("--atom-csv", type=Path, default=DEFAULT_ATOM_CSV)
    parser.add_argument("--out", type=Path, default=HERE / "results" / "sear_predictions.csv")
    args = parser.parse_args(argv)

    if not args.targets.exists():
        print(f"targets.json not found at {args.targets} — run prepare_opt.py first.", file=sys.stderr)
        return 1
    if not args.atom_csv.exists():
        print(f"atom_level.csv not found at {args.atom_csv} — run moldscript first.", file=sys.stderr)
        return 1

    targets = json.loads(args.targets.read_text())
    table = score(targets, args.atom_csv)
    if table.empty:
        print("No scored substrates — check that moldscript output covers the dataset.", file=sys.stderr)
        return 1

    args.out.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.out, index=False)

    print(f"\nScored {len(table)} substrates → {args.out.relative_to(HERE)}\n")
    summary_rows = []
    for desc in DESCRIPTORS:
        t1 = table[f"top1_{desc}"].dropna()
        t2 = table[f"top2_{desc}"].dropna()
        if t1.empty:
            continue
        summary_rows.append({
            "descriptor": desc,
            "top1_accuracy": round(t1.mean(), 3),
            "top2_accuracy": round(t2.mean(), 3),
            "n": int(t1.size),
        })
    summary = pd.DataFrame(summary_rows)
    print(summary.to_string(index=False))

    by_class = (
        table.groupby("class")[[f"top1_{d}" for d in DESCRIPTORS if f"top1_{d}" in table.columns]]
        .mean()
        .round(3)
    )
    print("\nTop-1 accuracy by substrate class:")
    print(by_class.to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
