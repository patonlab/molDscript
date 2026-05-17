"""Generate Fukui single-point inputs from completed opt logs.

Reads each `gaussian/opt/<name>_conf1.log`, extracts the final optimized
geometry via cclib, and writes three Gaussian SP inputs (neutral, +1, -1) at
that geometry with pop=npa so the moldscript fukui module can compute
condensed Fukui functions from natural charges.

File-naming matches moldscript's suffix convention (see arguments.txt):
  gaussian/neutral/<name>_conf1_neutral.com
  gaussian/oxidized/<name>_conf1_ox.com
  gaussian/reduced/<name>_conf1_red.com
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import cclib
import periodictable


HERE = Path(__file__).resolve().parent
DEFAULT_OPT_DIR = HERE / "gaussian" / "opt"
DEFAULT_OUT_ROOT = HERE / "gaussian"

SP_ROUTE = "#P B3LYP/6-31G(d) integral=ultrafine pop=npa"

CHARGE_STATES = [
    ("neutral", "neutral", 0, 1),
    ("oxidized", "ox", +1, 2),
    ("reduced", "red", -1, 2),
]


def extract_xyz(log_path: Path) -> list[tuple[str, float, float, float]]:
    data = cclib.io.ccread(str(log_path))
    if data is None or not hasattr(data, "atomcoords") or len(data.atomcoords) == 0:
        raise RuntimeError(f"cclib could not parse coordinates from {log_path}")
    coords = data.atomcoords[-1]
    nos = data.atomnos
    return [
        (str(periodictable.elements[int(z)]), float(x), float(y), float(z_))
        for z, (x, y, z_) in zip(nos, coords)
    ]


def write_sp(path: Path, xyz: list[tuple[str, float, float, float]], charge: int, mult: int, title: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write("%mem=4GB\n")
        f.write("%nprocshared=4\n")
        f.write(f"%chk={path.stem}.chk\n")
        f.write(f"{SP_ROUTE}\n\n")
        f.write(f"{title}\n\n")
        f.write(f"{charge} {mult}\n")
        for sym, x, y, z in xyz:
            f.write(f"{sym:>2}  {x:>12.6f}  {y:>12.6f}  {z:>12.6f}\n")
        f.write("\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--opt-dir", type=Path, default=DEFAULT_OPT_DIR)
    parser.add_argument("--out-root", type=Path, default=DEFAULT_OUT_ROOT)
    args = parser.parse_args(argv)

    logs = sorted(args.opt_dir.glob("*_conf1.log"))
    if not logs:
        print(f"No *_conf1.log found under {args.opt_dir} — run the opt jobs first.", file=sys.stderr)
        return 1

    for log in logs:
        name = log.stem.removesuffix("_conf1")
        xyz = extract_xyz(log)
        for sub_dir, suffix, charge, mult in CHARGE_STATES:
            out = args.out_root / sub_dir / f"{name}_conf1_{suffix}.com"
            write_sp(out, xyz, charge, mult, title=f"{name} {sub_dir} SP B3LYP/6-31G(d)")
        print(f"  {name}: wrote neutral / ox / red SPs")
    print(f"\nWrote SPs for {len(logs)} substrates under {args.out_root.relative_to(HERE)}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
