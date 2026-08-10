"""Compact steric descriptors for multi-frame XYZ ensembles.

The reader supports the concatenated XYZ files written by CREST and ORCA GOAT,
as well as any equivalent multi-frame XYZ supplied directly.  Per-conformer
steric values are calculated internally and summarized in molDscript's
existing molecule- and atom-level data dictionaries.
"""

from dataclasses import dataclass
import ast
import datetime
import math
from pathlib import Path
import re
import time
from typing import Tuple

import numpy as np
from dbstep.constants import bondi
from rdkit import Chem

from moldscript.argument_parser import load_variables
from moldscript.utils import molecule_keys, run_file_jobs


VDW_SCALE = 1.17
UNKNOWN_VDW_RADIUS = 2.0
MAX_GRID_CUBE_POINTS = 4_000_000
HARTREE_TO_KCAL_MOL = 627.5094740631
GAS_CONSTANT_KCAL_MOL_K = 0.00198720425864083

_FLOAT_PATTERN = (
    r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"
)
_LEADING_ENERGY_RE = re.compile(
    rf"^\s*({_FLOAT_PATTERN})(?:\s|$)", re.IGNORECASE
)
_LABELLED_ENERGY_RE = re.compile(
    rf"\benergy\s*(?:[:=]\s*|\s+)({_FLOAT_PATTERN})",
    re.IGNORECASE,
)
_UNCONVERGED_RE = re.compile(
    r"\bconverged\s*=\s*false\b", re.IGNORECASE
)


@dataclass(frozen=True)
class XYZEnsemble:
    """One conformer ensemble parsed from a concatenated XYZ file."""

    name: str
    path: Path
    symbols: Tuple[str, ...]
    atomic_numbers: np.ndarray
    coordinates: np.ndarray
    energies_hartree: np.ndarray

    @property
    def n_conformers(self):
        return int(self.coordinates.shape[0])

    @property
    def n_atoms(self):
        return int(self.coordinates.shape[1])


def _strip_user_suffix(name, suffix):
    suffix = str(suffix or "").strip().strip("_")
    if not suffix:
        return name
    token = f"_{suffix}"
    if name.endswith(token):
        return name[: -len(token)]
    if name == suffix:
        return ""
    return name


def _ensemble_name(path, suffix=""):
    path = Path(path)
    lower_name = path.name.lower()
    if lower_name == "crest_conformers.xyz":
        name = path.parent.name
    elif lower_name.endswith(".finalensemble.xyz"):
        name = path.name[: -len(".finalensemble.xyz")]
    else:
        name = path.stem
        marker = "_crest_conformers"
        if name.lower().endswith(marker):
            name = name[: -len(marker)]
    name = _strip_user_suffix(name, suffix)
    if not name:
        raise ValueError(
            f"Could not derive a molecule name from ensemble file {path}"
        )
    return name


def discover_ensemble_xyz(path):
    """Return supported ensemble XYZ files from a file or directory.

    Direct file input accepts any ``.xyz`` name.  Directory discovery is
    intentionally conservative to avoid treating trajectory, best-structure,
    or GOAT iteration files as independent ensembles.
    """

    path = Path(path).expanduser()
    if not path.exists():
        raise FileNotFoundError(f"Ensemble input does not exist: {path}")
    if path.is_file():
        if path.suffix.lower() != ".xyz":
            raise ValueError(
                f"Ensemble input must be an XYZ file, not {path.name}"
            )
        return [path.resolve()]

    selected = []
    for candidate in path.rglob("*.xyz"):
        name = candidate.name.lower()
        if (
            name == "crest_conformers.xyz"
            or name.endswith("_crest_conformers.xyz")
            or name.endswith(".finalensemble.xyz")
        ):
            selected.append(candidate.resolve())
    selected.sort(key=lambda item: str(item))
    if not selected:
        raise FileNotFoundError(
            "No standard ensemble XYZ files were found below "
            f"{path}. Expected crest_conformers.xyz, "
            "*_crest_conformers.xyz, or *.finalensemble.xyz; pass a "
            "nonstandard XYZ filename directly."
        )
    return selected


def _parse_comment_energy(comment, path, frame_number):
    if _UNCONVERGED_RE.search(comment):
        raise ValueError(
            f"{path}: frame {frame_number} is marked converged=false"
        )
    match = _LEADING_ENERGY_RE.search(comment)
    if match is None:
        match = _LABELLED_ENERGY_RE.search(comment)
    if match is None:
        raise ValueError(
            f"{path}: frame {frame_number} expected a Hartree energy on "
            "the XYZ comment line"
        )
    energy = float(match.group(1).replace("D", "E").replace("d", "e"))
    if not math.isfinite(energy):
        raise ValueError(
            f"{path}: frame {frame_number} ensemble energy must be finite"
        )
    return energy


def _normalise_symbol(value, path, line_number):
    symbol = value[:1].upper() + value[1:].lower()
    try:
        atomic_number = Chem.GetPeriodicTable().GetAtomicNumber(symbol)
    except RuntimeError as exc:
        raise ValueError(
            f"{path}: line {line_number} has unknown element symbol {value!r}"
        ) from exc
    if atomic_number <= 0:
        raise ValueError(
            f"{path}: line {line_number} has unknown element symbol {value!r}"
        )
    return symbol, atomic_number


def read_xyz_ensemble(path, suffix=""):
    """Parse a strict concatenated XYZ ensemble with Hartree energies."""

    path = Path(path).expanduser().resolve()
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    frame_coordinates = []
    energies = []
    expected_natoms = None
    expected_symbols = None
    atomic_numbers = None
    position = 0

    while position < len(lines):
        while position < len(lines) and not lines[position].strip():
            position += 1
        if position >= len(lines):
            break

        frame_number = len(energies) + 1
        count_line_number = position + 1
        try:
            natoms = int(lines[position].strip())
        except ValueError as exc:
            raise ValueError(
                f"{path}: line {count_line_number} expected an integer atom "
                "count"
            ) from exc
        if natoms <= 0:
            raise ValueError(
                f"{path}: frame {frame_number} atom count must be positive"
            )
        position += 1
        if position >= len(lines):
            raise ValueError(
                f"{path}: frame {frame_number} is missing its comment line"
            )

        energies.append(
            _parse_comment_energy(lines[position], path, frame_number)
        )
        position += 1

        symbols = []
        numbers = []
        coordinates = []
        for atom_index in range(natoms):
            if position >= len(lines):
                raise ValueError(
                    f"{path}: frame {frame_number} truncated after "
                    f"{atom_index} of {natoms} atom rows"
                )
            line_number = position + 1
            fields = lines[position].split()
            if len(fields) < 4:
                raise ValueError(
                    f"{path}: frame {frame_number} atom row {atom_index + 1} "
                    "requires an element and three coordinates"
                )
            symbol, atomic_number = _normalise_symbol(
                fields[0], path, line_number
            )
            try:
                xyz = [
                    float(value.replace("D", "E").replace("d", "e"))
                    for value in fields[1:4]
                ]
            except ValueError as exc:
                raise ValueError(
                    f"{path}: frame {frame_number} atom row {atom_index + 1} "
                    "has invalid coordinates"
                ) from exc
            if not np.all(np.isfinite(xyz)):
                raise ValueError(
                    f"{path}: frame {frame_number} coordinates must be finite"
                )
            symbols.append(symbol)
            numbers.append(atomic_number)
            coordinates.append(xyz)
            position += 1

        symbol_tuple = tuple(symbols)
        if expected_natoms is None:
            expected_natoms = natoms
            expected_symbols = symbol_tuple
            atomic_numbers = np.asarray(numbers, dtype=int)
        elif natoms != expected_natoms:
            raise ValueError(
                f"{path}: frame {frame_number} contains {natoms} atoms; "
                f"expected {expected_natoms}"
            )
        elif symbol_tuple != expected_symbols:
            raise ValueError(
                f"{path}: frame {frame_number} changes atom ordering or "
                "element identities"
            )
        frame_coordinates.append(coordinates)

    if not energies:
        raise ValueError(f"{path}: no XYZ conformers were found")

    return XYZEnsemble(
        name=_ensemble_name(path, suffix=suffix),
        path=path,
        symbols=expected_symbols,
        atomic_numbers=atomic_numbers,
        coordinates=np.asarray(frame_coordinates, dtype=float),
        energies_hartree=np.asarray(energies, dtype=float),
    )


def parse_radii(value):
    """Normalize a radius or list of radii to unique positive floats."""

    if value is None:
        value = [3.5]
    if isinstance(value, str):
        stripped = value.strip()
        if not stripped:
            return []
        try:
            value = ast.literal_eval(stripped)
        except (SyntaxError, ValueError):
            value = [part.strip() for part in stripped.split(",")]
    if np.isscalar(value):
        value = [value]

    radii = []
    for item in value:
        try:
            radius = float(item)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Ensemble radii must be numbers, not {item!r}"
            ) from exc
        if not math.isfinite(radius) or radius <= 0:
            raise ValueError("Ensemble radii must be finite and positive")
        if radius not in radii:
            radii.append(radius)
    return radii


def parse_atom_indices(value):
    """Parse 1-based atom indices and return their 0-based equivalents."""

    if value is None or value is False:
        return []
    if isinstance(value, str):
        stripped = value.strip()
        if not stripped:
            return []
        try:
            value = ast.literal_eval(stripped)
        except (SyntaxError, ValueError):
            value = [part.strip() for part in stripped.split(",")]
    if np.isscalar(value):
        value = [value]

    indices = []
    for item in value:
        try:
            numeric = float(item)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Ensemble atom indices must be whole numbers, not {item!r}"
            ) from exc
        if not math.isfinite(numeric) or not numeric.is_integer():
            raise ValueError("Ensemble atom indices must be whole numbers")
        index = int(numeric)
        if index < 1:
            raise ValueError("Ensemble atom indices are 1-based and positive")
        zero_based = index - 1
        if zero_based not in indices:
            indices.append(zero_based)
    return indices


def boltzmann_weights(energies_hartree, temperature=298.15):
    """Return normalized populations from conformer energies in Hartree."""

    energies = np.asarray(energies_hartree, dtype=float)
    if energies.ndim != 1 or energies.size == 0:
        raise ValueError("At least one ensemble energy is required")
    if not np.all(np.isfinite(energies)):
        raise ValueError("Ensemble energies must be finite")
    temperature = float(temperature)
    if not math.isfinite(temperature) or temperature <= 0:
        raise ValueError("Boltzmann temperature must be finite and positive")

    relative_kcal = (
        energies - np.min(energies)
    ) * HARTREE_TO_KCAL_MOL
    exponents = -relative_kcal / (
        GAS_CONSTANT_KCAL_MOL_K * temperature
    )
    exponents -= np.max(exponents)
    weights = np.exp(exponents)
    weights /= weights.sum()
    return weights


def _atomic_masses(atomic_numbers):
    table = Chem.GetPeriodicTable()
    return np.asarray(
        [table.GetAtomicWeight(int(number)) for number in atomic_numbers],
        dtype=float,
    )


def _atom_geometry_signature(conformer, atomic_numbers, atom_index):
    distances = np.linalg.norm(
        conformer - conformer[atom_index], axis=1
    )
    return tuple(
        sorted(
            (
                int(number),
                round(float(distance), 10),
            )
            for number, distance in zip(atomic_numbers, distances)
        )
    )


def _canonicalize_with_anchors(
    conformer, centered, atomic_numbers, tolerance
):
    """Fallback frame for linear or principal-axis-degenerate structures."""

    atom_indices = range(len(conformer))

    def first_key(atom_index):
        return (
            int(atomic_numbers[atom_index] != 1),
            round(float(np.linalg.norm(centered[atom_index])), 10),
            int(atomic_numbers[atom_index]),
            _atom_geometry_signature(
                conformer, atomic_numbers, atom_index
            ),
        )

    first_index = max(atom_indices, key=first_key)
    first_vector = centered[first_index]
    first_norm = np.linalg.norm(first_vector)
    if first_norm <= tolerance:
        return np.zeros_like(centered)
    first_vector = first_vector / first_norm

    def second_key(atom_index):
        cross_norm = np.linalg.norm(
            np.cross(first_vector, centered[atom_index])
        )
        return (
            round(float(cross_norm), 10),
            int(atomic_numbers[atom_index] != 1),
            int(atomic_numbers[atom_index]),
            _atom_geometry_signature(
                conformer, atomic_numbers, atom_index
            ),
        )

    second_index = max(atom_indices, key=second_key)
    second_vector = np.cross(
        first_vector, centered[second_index]
    )
    second_norm = np.linalg.norm(second_vector)
    if second_norm <= tolerance:
        result = np.zeros_like(centered)
        result[:, 0] = centered @ first_vector
        return result

    second_vector = second_vector / second_norm
    third_vector = np.cross(second_vector, first_vector)
    basis = np.column_stack(
        (first_vector, third_vector, second_vector)
    )
    return centered @ basis


def _canonicalize_coordinates(coordinates, atomic_numbers):
    """Place conformers in rotation- and atom-permutation-invariant frames."""

    canonical = np.zeros_like(coordinates, dtype=float)
    tolerance = 1.0e-10

    for conformer_index, conformer in enumerate(coordinates):
        centered = conformer - conformer.mean(axis=0)
        tensor = centered.T @ centered
        eigenvalues, basis = np.linalg.eigh(tensor)
        order = np.argsort(eigenvalues)[::-1]
        eigenvalues = eigenvalues[order]
        basis = basis[:, order]
        scale = max(float(eigenvalues[0]), tolerance)
        degenerate = np.any(
            np.abs(np.diff(eigenvalues)) <= scale * 1.0e-8
        )
        if degenerate:
            canonical[conformer_index] = _canonicalize_with_anchors(
                conformer,
                centered,
                atomic_numbers,
                tolerance,
            )
        else:
            canonical[conformer_index] = centered @ basis
    return canonical


def _shape_descriptors(coordinates, masses):
    radii_of_gyration = []
    anisotropies = []
    for conformer in coordinates:
        center_of_mass = np.sum(
            conformer * masses[:, None], axis=0
        ) / masses.sum()
        centered = conformer - center_of_mass
        tensor = (
            centered.T @ (centered * masses[:, None])
        ) / masses.sum()
        eigenvalues = np.clip(np.linalg.eigvalsh(tensor), 0.0, None)
        eigenvalue_sum = float(eigenvalues.sum())
        radii_of_gyration.append(math.sqrt(eigenvalue_sum))
        if eigenvalue_sum <= np.finfo(float).eps:
            anisotropies.append(0.0)
        else:
            mean = eigenvalue_sum / 3.0
            anisotropies.append(
                1.5
                * float(np.sum((eigenvalues - mean) ** 2))
                / (eigenvalue_sum**2)
            )
    return (
        np.asarray(radii_of_gyration),
        np.asarray(anisotropies),
    )


def _sphere_grid(max_radius, spacing):
    """Return an origin-anchored spherical voxel grid and squared radii."""

    max_radius = float(max_radius)
    spacing = float(spacing)
    if not math.isfinite(spacing) or spacing <= 0:
        raise ValueError("Ensemble grid spacing must be finite and positive")
    half_width = int(math.ceil(max_radius / spacing))
    side = 2 * half_width + 1
    cube_points = side**3
    if cube_points > MAX_GRID_CUBE_POINTS:
        raise ValueError(
            "Requested ensemble buried-volume grid exceeds the safety limit "
            f"of {MAX_GRID_CUBE_POINTS:,} cube points; increase "
            "--ensemble_grid or reduce --ensemble_radii"
        )
    axis = np.arange(-half_width, half_width + 1, dtype=float) * spacing
    grid = np.stack(
        np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1
    ).reshape(-1, 3)
    squared = np.einsum("ij,ij->i", grid, grid)
    inside = squared <= max_radius**2 + 1.0e-12
    return grid[inside], squared[inside]


def _vdw_radii(symbols):
    missing = sorted({symbol for symbol in symbols if symbol not in bondi})
    radii = np.asarray(
        [bondi.get(symbol, UNKNOWN_VDW_RADIUS) * VDW_SCALE for symbol in symbols],
        dtype=float,
    )
    return radii, missing


def _buried_volume_series(
    coordinates,
    centers,
    symbols,
    center_index,
    radii,
    grid_spacing,
    include_hydrogen,
    excluded_indices,
    grid_data=None,
):
    atom_count = len(symbols)
    invalid = [index + 1 for index in excluded_indices if index >= atom_count]
    if invalid:
        raise ValueError(
            f"Ensemble buried-volume exclusions {invalid} are outside the "
            f"1-{atom_count} atom range"
        )

    included = np.ones(atom_count, dtype=bool)
    if not include_hydrogen:
        included &= np.asarray(symbols) != "H"
    included[center_index] = False
    if excluded_indices:
        included[np.asarray(excluded_indices, dtype=int)] = False

    included_symbols = np.asarray(symbols)[included]
    selected_radii, missing = _vdw_radii(included_symbols)
    if grid_data is None:
        grid_data = _sphere_grid(max(radii), grid_spacing)
    points, squared_from_center = grid_data
    sphere_masks = {
        radius: squared_from_center <= radius**2 + 1.0e-12
        for radius in radii
    }
    values = {
        radius: np.zeros(len(coordinates), dtype=float) for radius in radii
    }

    for conformer_index, (conformer, center_point) in enumerate(
        zip(coordinates, centers)
    ):
        relative_atoms = conformer[included] - center_point
        occupied = np.zeros(len(points), dtype=bool)
        for atom_coordinate, atom_radius in zip(
            relative_atoms, selected_radii
        ):
            delta = points - atom_coordinate
            occupied |= (
                np.einsum("ij,ij->i", delta, delta)
                <= atom_radius**2 + 1.0e-12
            )
        for radius, sphere_mask in sphere_masks.items():
            denominator = int(np.count_nonzero(sphere_mask))
            values[radius][conformer_index] = (
                100.0
                * np.count_nonzero(occupied & sphere_mask)
                / denominator
            )
    return values, missing


def _buried_volume_by_atom(
    coordinates,
    symbols,
    radii,
    grid_spacing,
    include_hydrogen,
    excluded_indices,
):
    """Return per-conformer buried volume around every atom."""

    atom_count = len(symbols)
    grid_data = _sphere_grid(max(radii), grid_spacing)
    volume_matrices = {
        radius: np.zeros(
            (len(coordinates), atom_count), dtype=float
        )
        for radius in radii
    }
    missing_radii = set()

    for center_index in range(atom_count):
        center_values, missing = _buried_volume_series(
            coordinates,
            coordinates[:, center_index, :],
            symbols,
            center_index,
            radii,
            grid_spacing,
            include_hydrogen,
            excluded_indices,
            grid_data=grid_data,
        )
        missing_radii.update(missing)
        for radius, values in center_values.items():
            volume_matrices[radius][:, center_index] = values

    return volume_matrices, sorted(missing_radii)


def _add_min_max_range(target, root, values, unit=None):
    values = np.asarray(values, dtype=float)
    suffix = f"_{unit}" if unit else ""
    minimum = float(np.min(values))
    maximum = float(np.max(values))
    target[f"{root}_min{suffix}"] = minimum
    target[f"{root}_max{suffix}"] = maximum
    target[f"{root}_range{suffix}"] = maximum - minimum


def _radius_label(radius):
    return f"{float(radius):g}".replace(".", "_")


def _analyze_ensemble_job(job):
    (
        source_path,
        suffix,
        radii,
        grid_spacing,
        temperature,
        include_hydrogen,
        excluded_indices,
    ) = job
    ensemble_data = read_xyz_ensemble(source_path, suffix=suffix)
    lowest_index = int(np.argmin(ensemble_data.energies_hartree))
    canonical_coordinates = _canonicalize_coordinates(
        ensemble_data.coordinates,
        ensemble_data.atomic_numbers,
    )
    masses = _atomic_masses(ensemble_data.atomic_numbers)

    radius_of_gyration, shape_anisotropy = _shape_descriptors(
        canonical_coordinates, masses
    )
    mol_descriptors = {}
    _add_min_max_range(
        mol_descriptors,
        "ensemble_radius_of_gyration",
        radius_of_gyration,
        unit="angstrom",
    )
    _add_min_max_range(
        mol_descriptors,
        "ensemble_shape_anisotropy",
        shape_anisotropy,
    )

    atom_descriptors = {}
    missing_radii = []
    if radii:
        weights = boltzmann_weights(
            ensemble_data.energies_hartree,
            temperature=temperature,
        )
        volume_matrices, missing_radii = _buried_volume_by_atom(
            canonical_coordinates,
            ensemble_data.symbols,
            radii,
            grid_spacing,
            include_hydrogen,
            excluded_indices,
        )
        for radius, values in volume_matrices.items():
            root = (
                f"ensemble_buried_volume_r_"
                f"{_radius_label(radius)}A"
            )
            atom_descriptors[f"{root}_min_percent"] = np.min(
                values, axis=0
            )
            atom_descriptors[f"{root}_max_percent"] = np.max(
                values, axis=0
            )
            atom_descriptors[
                f"{root}_boltzmann_mean_percent"
            ] = np.sum(values * weights[:, None], axis=0)
            atom_descriptors[
                f"{root}_lowest_energy_percent"
            ] = values[lowest_index].copy()

    return {
        "name": ensemble_data.name,
        "path": ensemble_data.path,
        "atomic_numbers": ensemble_data.atomic_numbers,
        "minimum_energy_hartree": float(
            ensemble_data.energies_hartree[lowest_index]
        ),
        "mol": mol_descriptors,
        "atom": atom_descriptors,
        "missing_radii": missing_radii,
    }


def _merge_ensemble_result(data_dict, result, require_existing=False):
    name = result["name"]
    available = molecule_keys(data_dict)
    if require_existing and name not in available:
        raise ValueError(
            f"Could not match ensemble key {name!r} to an existing molecule. "
            f"Existing molecule keys: {', '.join(available) or 'none'}. "
            "Pass --suffix_ensemble when the ensemble filename contains an "
            "extra trailing tag."
        )

    if name in available:
        entry = data_dict[name]
        existing_numbers = np.asarray(
            entry["atom"].get("atomnos"), dtype=int
        )
        if not np.array_equal(
            existing_numbers, result["atomic_numbers"]
        ):
            raise ValueError(
                f"Ensemble {name!r} does not have the same atom count and "
                "element order as the existing molecular data"
            )
    else:
        entry = {
            "mol": {
                "smiles": "",
                "scfenergy": result["minimum_energy_hartree"],
            },
            "atom": {
                "atomnos": result["atomic_numbers"],
            },
            "bond": {},
            "CPU_time": datetime.timedelta(0),
        }
        data_dict[name] = entry

    entry["mol"].update(result["mol"])
    entry["atom"].update(result["atom"])
    entry.setdefault("CPU_time", datetime.timedelta(0))


class ensemble:
    """Analyze generic XYZ ensembles and update a molDscript data dictionary."""

    def __init__(
        self,
        path,
        data_dict=None,
        radii=None,
        grid=0.25,
        temp=298.15,
        include_h=False,
        exclude="",
        suffix="",
        output="",
        workers=1,
        create_dat=True,
        **kwargs,
    ):
        started = time.time()
        options = {
            "output": output,
            "workers": workers,
        }
        options.update(kwargs)
        self.args = load_variables(
            options, "ENSEMBLE", create_dat=create_dat
        )
        self.path = path
        self.data_dict = data_dict if data_dict is not None else {}
        self.file_data = self.data_dict
        self.radii = parse_radii(radii)
        self.excluded_indices = parse_atom_indices(exclude)
        self.ensemble_files = discover_ensemble_xyz(path)

        names = [
            _ensemble_name(item, suffix=suffix)
            for item in self.ensemble_files
        ]
        duplicates = sorted(
            {name for name in names if names.count(name) > 1}
        )
        if duplicates:
            raise ValueError(
                "Multiple ensemble XYZ files map to the same molecule key: "
                + ", ".join(duplicates)
            )

        self.args.log.write(
            "-- XYZ Ensemble Parameter Collection starting"
        )
        original_molecule_keys = set(molecule_keys(self.data_dict))
        self.data_dict.setdefault("CPU_time", [])
        jobs = [
            (
                source_path,
                suffix,
                self.radii,
                float(grid),
                float(temp),
                bool(include_h),
                self.excluded_indices,
            )
            for source_path in self.ensemble_files
        ]
        self.results = run_file_jobs(
            jobs,
            _analyze_ensemble_job,
            workers=workers,
            logger=self.args.log,
        )

        for result in self.results:
            _merge_ensemble_result(
                self.data_dict,
                result,
                require_existing=bool(original_molecule_keys),
            )
            details = (
                f"o  Added compact ensemble descriptors from "
                f"{result['path'].name} to {result['name']}"
            )
            self.args.log.write_only(details)
            if self.radii:
                atom_count = len(result["atomic_numbers"])
                self.args.log.write_only(
                    f"   Buried volume evaluated around all "
                    f"{atom_count} atoms"
                )
            if result["missing_radii"]:
                self.args.log.write(
                    "Warning: DBSTEP Bondi radii were unavailable for "
                    f"{', '.join(result['missing_radii'])}; used "
                    f"{UNKNOWN_VDW_RADIUS:.1f} Angstrom before the "
                    f"{VDW_SCALE:.2f} scale factor"
                )

        elapsed = round(time.time() - started, 2)
        self.args.log.write(
            "-- XYZ Ensemble Parameter Collection complete in "
            f"{elapsed} seconds"
        )
        if create_dat:
            self.args.log.finalize()


Ensemble = ensemble
