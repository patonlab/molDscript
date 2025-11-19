######################################################.
#          This file stores functions used           #
#                in multiple modules                 #
######################################################.

import os
import ast
from pathlib import Path
import glob
import datetime
import numpy as np
import cclib as cc
import moldscript.xyz2mol as xyz2mol
from rdkit import Chem

k_B_hartree = 3.1668114e-6  # hartree/K
J_TO_AU = 4.184 * 627.509541 * 1000.0  # UNIT CONVERSION
eV_to_hartree = 0.0367493

# class for logging
class Logger:
    """Simple file logger used to emit MOLDSCRIPT_*.dat audit files."""

    def __init__(self, file_path=None, verbose=True):
        self.verbose = verbose
        self.path = Path(file_path) if file_path else None
        self._handle = None
        if self.verbose and self.path:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            self._handle = self.path.open("w", encoding="utf-8")

    @classmethod
    def silent(cls):
        """Return a logger that swallows all output."""
        return cls(file_path=None, verbose=False)

    def write(self, message):
        """Write a message to the log file (if enabled) and echo to stdout."""
        if self._handle:
            self._handle.write(f"{message}\n")
            self._handle.flush()
        print(message)

    def write_only(self, message):
        """Write a message only to the log file (if enabled)."""
        if self._handle:
            self._handle.write(f"{message}\n")
            self._handle.flush()

    def finalize(self):
        """Close the underlying file handle, if one was opened."""
        if self._handle:
            self._handle.close()
            self._handle = None


def build_log_path(output_prefix: str, module_code: str, suffix: str = "dat") -> Path:
    """Return the filesystem path for a module audit log."""
    prefix = output_prefix or ""
    filename = f"{prefix}MOLDSCRIPT_{module_code}.{suffix}"
    path = Path(filename)
    if not path.is_absolute():
        path = Path.cwd() / path
    return path


def _normalise_cpu_spans(cpu_times):
    if cpu_times is None:
        return []
    if isinstance(cpu_times, (datetime.timedelta, int, float)):
        cpu_times = [cpu_times]
    spans = []
    for span in cpu_times:
        if span is None:
            continue
        if isinstance(span, datetime.timedelta):
            spans.append(span)
        elif isinstance(span, (int, float)):
            spans.append(datetime.timedelta(seconds=float(span)))
        else:
            try:
                spans.append(datetime.timedelta(seconds=float(span)))
            except (TypeError, ValueError):
                continue
    return spans


def record_cpu_time(data_dict, file_key, source_path, cpu_times):
    """Accumulate CPU times for a parsed file and avoid double counting."""
    spans = _normalise_cpu_spans(cpu_times)
    if not spans:
        return 0.0

    bucket = None
    if source_path is not None:
        bucket = data_dict.setdefault("CPU_time", [])
        if source_path in bucket:
            return 0.0

    entry = data_dict.setdefault(file_key, {})
    total = entry.get("CPU_time", datetime.timedelta(0))
    added = datetime.timedelta(0)
    for delta in spans:
        total += delta
        added += delta
    entry["CPU_time"] = total

    if source_path is not None and bucket is not None:
        bucket.append(source_path)

    return added.total_seconds()


def format_timedelta(value: datetime.timedelta) -> str:
    total_seconds = value.total_seconds()
    total_hours = total_seconds / 3600.0
    return f"{total_hours:.2f} hours"


def add_cpu_times(file_data):
    ''' add cpu times for all files'''
    total_cpu = datetime.timedelta(0)
    for entry in file_data.values():
        if not isinstance(entry, dict):
            continue
        cpu_value = entry.get("CPU_time") or entry.get("cpu_time")
        if cpu_value:
            total_cpu += cpu_value
    return total_cpu
def initiate_data_dict(data, logger=None):
    """
    Initiates a data dictionary to store all the data from the files.
    Progress and informational messages are emitted to `logger` when provided,
    otherwise they fall back to printing to stdout.
    """
    if logger:
        logger.write(f"Initializing data parsing with SMILES and geometry data")
    else:
        print(f"Initializing data parsing with SMILES and geometry data")

    total = len(data)
    data_dict = {}
    data_dict["CPU_time"] = []

    if total == 0:
        if logger:
            logger.write("No files to process.")
        else:
            print("No files to process.")
        return data_dict

    last_step = 0  # tracks 5% steps printed (0..20)
    for i, file_name in enumerate(data.keys(), start=1):
        # Progress reporting every 5%
        percent = int((i / total) * 100)
        step = percent // 5
        if step > last_step:
            # emit any missed intermediate 5% markers if step jumped
            for s in range(last_step + 1, step + 1):
                if logger:
                    logger.write(f"Progress: {s * 5}% ({i}/{total})")
                else:
                    print(f"Progress: {s * 5}% ({i}/{total})")
            last_step = step

        data_dict[file_name] = dict()
        data_dict[file_name]["mol"] = dict()
        data_dict[file_name]["atom"] = dict()
        data_dict[file_name]["bond"] = dict()
        parsed_data = parse_cc_data(file_name, data[file_name])
        try:
            mol = xyz2mol.xyz2mol(parsed_data.atomnos.tolist(), parsed_data.atomcoords[-1].tolist(), charge=parsed_data.charge)[0]
            smi = Chem.MolToSmiles(mol)
        except:
            if logger:
                logger.write("Encountered an issue with the mol embedding. Skipping smiles string.")
            else:
                print("Encountered an issue with the mol embedding. Skipping smiles string.")
            smi = ''
        data_dict[file_name]["mol"]["smiles"] = smi
        data_dict[file_name]["atom"]["atomnos"] = parsed_data.atomnos
        data_dict[file_name]["bond"]["bond_length"] = parsed_data.bond_data_matrix
        data_dict[file_name]["mol"]["scfenergy"] = (parsed_data.scfenergies[-1] * eV_to_hartree)

    # Ensure 100% is emitted
    if last_step < 20:
        for s in range(last_step + 1, 21):
            if logger:
                logger.write(f"Progress: {s * 5}% ({total}/{total})")
            else:
                print(f"Progress: {s * 5}% ({total}/{total})")

    return data_dict

def format_lists(value):
    '''
    Transforms strings into a list
    '''

    if not isinstance(value, list):
        try:
            value = ast.literal_eval(value)
        except (SyntaxError, ValueError):
            # this line fixes issues when using "[X]" or ["X"] instead of "['X']" when using lists
            value = value.replace('[',']').replace(',',']').replace("'",']').split(']')
            while('' in value):
                value.remove('')
    return value
def bond_data_matrix(data):
        try:
            coords = data.atomcoords[-1]
        except:
            coords = data
        bond_data_matrix_list = []
        for atom1 in range(len(coords)):
            row = []
            for atom2 in range(len(coords)):
                p1 = np.array(coords[atom1])
                p2 = np.array(coords[atom2])
                squared_dist = np.sum((p1 - p2) ** 2, axis=0)
                dist = np.sqrt(squared_dist)
                row.append(dist)
            bond_data_matrix_list.append(row)
        return bond_data_matrix_list
def parse_cc_data(file_name, file):
        try:
            parser = cc.io.ccopen(file)
            cc_data = parser.parse()
            setattr(cc_data, "bond_data_matrix", bond_data_matrix(cc_data))
        except:
            print(
                f"\nx  Could not parse {file_name}")
            raise SystemExit(f"Error parsing {file_name}. Ensure the file is a valid cclib file.")
        return cc_data
def get_files(value):

        if value[-1]=='/':
            value = value[:-1]
        if (Path(f"{value}").exists() and os.getcwd() not in f"{value}"):
            list_of_val_log = glob.glob(f"{os.getcwd()}/{value}/*.log")
            list_of_val_out = glob.glob(f"{os.getcwd()}/{value}/*.out")
        else:
            list_of_val_log = glob.glob(f"{value}/*.log")
            list_of_val_out = glob.glob(f"{value}/*.out")
        length_out = len(list_of_val_out)
        length_log = len(list_of_val_log)
        if length_log >= length_out:
            list_of_val = list_of_val_log
        else:
            list_of_val = list_of_val_out
        return list_of_val

def find_nth(haystack: str, needle: str, n: int) -> int:
    start = haystack.find(needle)
    while start >= 0 and n > 1:
        start = haystack.find(needle, start+len(needle))
        n -= 1
    return start
def get_filename(fullname, dd):
    flist = list(dd.keys())
    tempname = fullname
    for i in range(fullname.count("_")+1):
        try:
            findex = flist.index(tempname)
            keyname = flist[findex]
            return keyname
        except:
            tempname = tempname.rsplit("_", 1)[0]
            print(tempname)
    print(
        f"Error processing file {fullname}. Ensure consistent naming as described in the docs."
    )
    raise SystemExit


