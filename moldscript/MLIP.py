######################################################.
#        This file stores the mlip class             #
######################################################.

import datetime
import shlex
from pathlib import Path

import numpy as np
from rdkit import Chem

from moldscript.utils import eV_to_hartree


class mlip:
    """
    Parse MACE-Polar extxyz files into a dictionary compatible with get_df.

    Expected folder layout:
      - neutral/*.extxyz
      - reduced/*.extxyz (optional)
      - oxidized/*.extxyz (optional)

    The generated dictionary follows the same structure used by the other
    modules: {
      "CPU_time": [],
      "<name>": {
          "mol": {...},
          "atom": {...},
          "bond": {...},
          "CPU_time": datetime.timedelta(...)
      }
    }
    """

    def __init__(
        self,
        neutral,
        reduced="",
        oxidized="",
        data_dict=None,
        atom_charge_property="q_mace_polar",
        file_glob="*.extxyz",
    ):
        self.neutral = Path(neutral)
        self.reduced = Path(reduced) if reduced else None
        self.oxidized = Path(oxidized) if oxidized else None
        self.file_glob = file_glob
        self.atom_charge_property = atom_charge_property
        self.ptable = Chem.GetPeriodicTable()

        self._progress_total = self._count_input_files()
        self._progress_done = 0
        self._progress_step = 0

        print("-- MLIP Parameter Collection starting", flush=True)

        if data_dict is None or data_dict == {}:
            self.data_dict = {"CPU_time": []}
        else:
            self.data_dict = data_dict
            self.data_dict.setdefault("CPU_time", [])

        # key -> {"neutral": rec, "reduced": rec, "oxidized": rec}
        self.state_records = {}

        self._parse_state_dir(self.neutral, "neutral")
        if self.reduced is not None:
            self._parse_state_dir(self.reduced, "reduced")
        if self.oxidized is not None:
            self._parse_state_dir(self.oxidized, "oxidized")

        self.file_data = self._build_data_dict()
        print("-- MLIP Parameter Collection complete", flush=True)

    def _count_input_files(self):
        total = len(list(self.neutral.glob(self.file_glob)))
        if self.reduced is not None:
            total += len(list(self.reduced.glob(self.file_glob)))
        if self.oxidized is not None:
            total += len(list(self.oxidized.glob(self.file_glob)))
        return total

    def _report_progress(self):
        if self._progress_total == 0:
            return
        percent = int((self._progress_done / self._progress_total) * 100)
        step = percent // 5
        if step > self._progress_step:
            for marker in range(self._progress_step + 1, step + 1):
                print(
                    f"Progress: {marker * 5}% ({self._progress_done}/{self._progress_total})",
                    flush=True,
                )
            self._progress_step = step

    def _parse_state_dir(self, state_dir, state):
        files = sorted(state_dir.glob(self.file_glob))
        total = len(files)
        for file_path in files:
            record = self._read_extxyz(file_path)
            key = self._canonical_key(record["header"], file_path)
            if key not in self.state_records:
                self.state_records[key] = {}
            self.state_records[key][state] = record
            self._progress_done += 1
            self._report_progress()
        if total == 0:
            print(f"!  No MLIP files found in {state_dir}", flush=True)

    def _build_data_dict(self):
        total = len(self.state_records)
        if total == 0:
            print("x  Could not find MLIP files to obtain information", flush=True)
            return self.data_dict

        for key, states in self.state_records.items():
            base_state = self._pick_base_state(states)
            base_record = states[base_state]

            symbols = base_record["symbols"]
            coords = np.array(base_record["coords"], dtype=float)
            atomnos = np.array([self.ptable.GetAtomicNumber(sym) for sym in symbols], dtype=int)

            entry = {
                "mol": {},
                "atom": {},
                "bond": {},
                "CPU_time": datetime.timedelta(0),
            }

            entry["atom"]["atomnos"] = atomnos
            entry["bond"]["bond_length"] = self._bond_length_matrix(coords)

            entry["mol"].setdefault("smiles", "")

            # State energies
            e_neutral_ev = self._state_energy_ev(states.get("neutral"))
            e_reduced_ev = self._state_energy_ev(states.get("reduced"))
            e_oxidized_ev = self._state_energy_ev(states.get("oxidized"))

            if np.isfinite(e_neutral_ev):
                entry["mol"]["scfenergy"] = e_neutral_ev * eV_to_hartree
            elif np.isfinite(e_reduced_ev):
                entry["mol"]["scfenergy"] = e_reduced_ev * eV_to_hartree
            elif np.isfinite(e_oxidized_ev):
                entry["mol"]["scfenergy"] = e_oxidized_ev * eV_to_hartree
            else:
                entry["mol"]["scfenergy"] = np.nan

            if np.isfinite(e_oxidized_ev) and np.isfinite(e_neutral_ev):
                ie_h = (e_oxidized_ev - e_neutral_ev) * eV_to_hartree
                entry["mol"]["vertical_ie"] = ie_h
                entry["mol"]["ionization_potential"] = ie_h
            else:
                entry["mol"]["vertical_ie"] = np.nan
                entry["mol"]["ionization_potential"] = np.nan

            if np.isfinite(e_reduced_ev) and np.isfinite(e_neutral_ev):
                ea_h = (e_reduced_ev - e_neutral_ev) * eV_to_hartree
                entry["mol"]["vertical_ea"] = ea_h
                entry["mol"]["electron_affinity"] = ea_h
            else:
                entry["mol"]["vertical_ea"] = np.nan
                entry["mol"]["electron_affinity"] = np.nan

            # Dipole from the first state that provides it (neutral preferred).
            dip = self._dipole_vector(states)
            entry["mol"]["dipole_x"] = dip[0]
            entry["mol"]["dipole_y"] = dip[1]
            entry["mol"]["dipole_z"] = dip[2]
            entry["mol"]["dipole"] = float(np.linalg.norm(dip)) if np.all(np.isfinite(dip)) else np.nan

            q_neutral = self._state_charges(states.get("neutral"), len(atomnos))
            q_reduced = self._state_charges(states.get("reduced"), len(atomnos))
            q_oxidized = self._state_charges(states.get("oxidized"), len(atomnos))

            # Follow sign convention already used in fukui.py
            if np.all(np.isfinite(q_neutral)) and np.all(np.isfinite(q_reduced)):
                fplus = -1.0 * (q_reduced - q_neutral)
            else:
                fplus = np.full(len(atomnos), np.nan)

            if np.all(np.isfinite(q_neutral)) and np.all(np.isfinite(q_oxidized)):
                fminus = -1.0 * (q_neutral - q_oxidized)
            else:
                fminus = np.full(len(atomnos), np.nan)

            frad = 0.5 * (fplus + fminus)
            entry["atom"]["fplus"] = fplus
            entry["atom"]["fminus"] = fminus
            entry["atom"]["frad"] = frad
            entry["atom"]["charges_neutral"] = q_neutral

            self.data_dict[key] = entry

        return self.data_dict

    def _read_extxyz(self, file_path):
        lines = file_path.read_text(encoding="utf-8").splitlines()
        if len(lines) < 2:
            raise ValueError(f"Invalid extxyz file: {file_path}")

        natoms = int(lines[0].strip())
        header_line = lines[1].strip()
        header = self._parse_header_line(header_line)

        schema = self._parse_properties_schema(header.get("Properties", "species:S:1:pos:R:3"))
        atom_lines = lines[2 : 2 + natoms]

        symbols = []
        coords = []
        atom_props = {name: [] for (name, ncols) in schema if name not in ("species", "pos") and ncols == 1}

        for line in atom_lines:
            toks = line.split()
            idx = 0
            row = {}
            for name, ncols in schema:
                vals = toks[idx : idx + ncols]
                idx += ncols
                row[name] = vals if ncols > 1 else vals[0]

            symbols.append(str(row["species"]))
            coords.append([float(x) for x in row["pos"]])

            for prop_name in atom_props.keys():
                atom_props[prop_name].append(self._to_scalar(row[prop_name]))

        for prop_name, values in atom_props.items():
            atom_props[prop_name] = np.array(values, dtype=float)

        return {
            "file_path": str(file_path),
            "natoms": natoms,
            "header": header,
            "symbols": symbols,
            "coords": coords,
            "atom_props": atom_props,
        }

    def _parse_header_line(self, line):
        out = {}
        for token in shlex.split(line):
            if "=" not in token:
                continue
            k, v = token.split("=", 1)
            out[k] = v
        return out

    def _parse_properties_schema(self, properties):
        # e.g. species:S:1:pos:R:3:q_mace_polar:R:1
        parts = properties.split(":")
        schema = []
        i = 0
        while i + 2 < len(parts):
            name = parts[i]
            ncols = int(parts[i + 2])
            schema.append((name, ncols))
            i += 3
        return schema

    def _canonical_key(self, header, file_path):
        source = header.get("source_file", "")
        if source:
            return Path(source).stem

        stem = file_path.stem
        suffixes = [
            "_neutral_vertical",
            "_anion_Nplus1_vertical",
            "_cation_Nminus1_vertical",
            "_neutral",
            "_anion_Nplus1",
            "_cation_Nminus1",
        ]
        for suffix in suffixes:
            if stem.endswith(suffix):
                return stem[: -len(suffix)]
        return stem

    def _pick_base_state(self, states):
        if "neutral" in states:
            return "neutral"
        if "reduced" in states:
            return "reduced"
        return "oxidized"

    def _state_energy_ev(self, record):
        if record is None:
            return np.nan
        return self._to_scalar(record["header"].get("mace_polar_energy_eV", np.nan))

    def _state_charges(self, record, natoms):
        if record is None:
            return np.full(natoms, np.nan)
        values = record["atom_props"].get(self.atom_charge_property)
        if values is None:
            return np.full(natoms, np.nan)
        arr = np.array(values, dtype=float)
        if len(arr) == natoms:
            return arr
        padded = np.full(natoms, np.nan)
        ncopy = min(natoms, len(arr))
        padded[:ncopy] = arr[:ncopy]
        return padded

    def _dipole_vector(self, states):
        for state in ["neutral", "reduced", "oxidized"]:
            record = states.get(state)
            if record is None:
                continue
            header = record["header"]
            x = self._to_scalar(header.get("mace_polar_dipole_x", np.nan))
            y = self._to_scalar(header.get("mace_polar_dipole_y", np.nan))
            z = self._to_scalar(header.get("mace_polar_dipole_z", np.nan))
            vec = np.array([x, y, z], dtype=float)
            if np.any(np.isfinite(vec)):
                return vec
        return np.array([np.nan, np.nan, np.nan], dtype=float)

    def _bond_length_matrix(self, coords):
        nat = len(coords)
        out = np.zeros((nat, nat), dtype=float)
        for i in range(nat):
            delta = coords - coords[i]
            out[i, :] = np.sqrt(np.sum(delta * delta, axis=1))
        return out

    def _to_scalar(self, value):
        try:
            return float(value)
        except (TypeError, ValueError):
            return value
