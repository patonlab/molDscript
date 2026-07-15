######################################################.
#          This file stores functions used           #
#                in multiple modules                 #
######################################################.

import os
import ast
import re
from pathlib import Path
import glob
import datetime
import numpy as np
import cclib as cc
import moldscript.xyz2mol as xyz2mol
from rdkit import Chem
from concurrent.futures import ProcessPoolExecutor, as_completed

try:
    from rich.console import Console
    from rich.panel import Panel
    from rich.progress import (
        BarColumn,
        MofNCompleteColumn,
        Progress,
        TextColumn,
        TimeElapsedColumn,
        TimeRemainingColumn,
    )
except ImportError:  # pragma: no cover - fallback for editable checkouts before deps install
    Console = None
    Panel = None
    Progress = None

try:
    from tqdm.auto import tqdm
except ImportError:  # pragma: no cover - final fallback for minimal environments
    tqdm = None

k_B_hartree = 3.1668114e-6  # hartree/K
J_TO_AU = 4.184 * 627.509541 * 1000.0  # UNIT CONVERSION
eV_to_hartree = 0.0367493

_LOG_PATHS_INITIALIZED = set()
_RUN_LOG_PATH = None
_PROGRESS_RE = re.compile(r"Progress:\s*\d+%\s*\((\d+)\s*/\s*(\d+)\)")


class Terminal:
    """Terminal output backed by Rich, then tqdm, then plain text."""

    def __init__(self):
        self.console = Console(highlight=False) if Console else None
        self.progress = None
        self.progress_backend = None
        self.task_id = None
        self.task_description = "Processing files"
        self.pending_description = "Processing files"

    def print(self, message="", style=None):
        if self.console:
            self.console.print(message, style=style)
        elif tqdm and self.progress_backend == "tqdm" and self.progress is not None:
            tqdm.write(str(message))
        else:
            print(message)

    def panel(self, message, title=None, border_style="cyan"):
        if self.console and Panel:
            self.console.print(
                Panel.fit(message, title=title, border_style=border_style)
            )
        else:
            if title:
                self.print(title)
            self.print(message)

    def set_pending_progress(self, description):
        if description:
            self.pending_description = description

    def update_progress(self, current, total):
        if total <= 0:
            return

        description = self.pending_description or "Processing files"
        if Progress:
            if self.progress is None or description != self.task_description:
                self.finish_progress()
                self.task_description = description
                self.progress_backend = "rich"
                self.progress = Progress(
                    TextColumn("[bold cyan]{task.description}"),
                    BarColumn(),
                    MofNCompleteColumn(),
                    TextColumn("{task.percentage:>3.0f}%"),
                    TimeElapsedColumn(),
                    TimeRemainingColumn(),
                    console=self.console,
                    transient=False,
                )
                self.progress.start()
                self.task_id = self.progress.add_task(description, total=total)

            self.progress.update(self.task_id, total=total, completed=current)
        elif tqdm:
            if self.progress is None or description != self.task_description:
                self.finish_progress()
                self.task_description = description
                self.progress_backend = "tqdm"
                self.progress = tqdm(total=total, desc=description, unit="file")

            if self.progress.total != total:
                self.progress.total = total
            delta = current - self.progress.n
            if delta > 0:
                self.progress.update(delta)
            elif delta < 0:
                self.progress.n = current
                self.progress.refresh()
        else:
            return

        if current >= total:
            self.finish_progress()

    def finish_progress(self):
        if self.progress is not None:
            if self.progress_backend == "rich":
                self.progress.stop()
            elif self.progress_backend == "tqdm":
                self.progress.close()
            self.progress = None
            self.progress_backend = None
            self.task_id = None


terminal = Terminal()


def _log_key(path):
    try:
        return str(Path(path).resolve())
    except OSError:
        return str(Path(path).absolute())


def append_run_log(message):
    """Append a line to the active run log, if one has been configured."""
    if _RUN_LOG_PATH is None:
        return
    path = Path(_RUN_LOG_PATH)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(f"{message}\n")


def terminal_print(message="", style=None):
    terminal.print(message, style=style)


def terminal_info(message):
    terminal.print(message, style="cyan")


def terminal_success(message):
    terminal.finish_progress()
    terminal.print(message, style="bold green")


def terminal_warning(message):
    terminal.print(message, style="yellow")


def terminal_error(message):
    terminal.finish_progress()
    terminal.print(message, style="bold red")


def emit(message, style=None, log=True):
    """Write a concise status line to the terminal and the active run log."""
    if log:
        append_run_log(message)
    terminal.print(message, style=style)


def print_run_header(version, timestamp, reference, argv=None):
    args = " ".join(argv or [])
    body = f"{timestamp}\n{reference}"
    if args:
        body += f"\n\nArguments: {args}"
    terminal.panel(body, title=f"molDscript v {version}", border_style="cyan")


def _display_log_message(message):
    text = str(message).strip()
    if not text:
        return

    progress_match = _PROGRESS_RE.search(text)
    if progress_match:
        current, total = [int(value) for value in progress_match.groups()]
        terminal.update_progress(current, total)
        return

    lowered = text.lower()
    if lowered.startswith("command line used"):
        return
    if lowered.startswith(("functional used", "basis set used", "package used")):
        return
    if " version used" in lowered or lowered.startswith("charges used"):
        return
    if lowered.startswith("initializing data parsing"):
        description = "Building baseline structures and geometry descriptors"
        terminal.set_pending_progress(description)
        return

    if "starting" in lowered:
        description = _clean_task_description(text)
        terminal.set_pending_progress(description)
        return

    if "complete" in lowered or "finished" in lowered:
        terminal_success(_clean_completion_message(text))
        return

    if text.startswith("x") or "error" in lowered:
        terminal_error(text)
    elif text.startswith("!") or "warning" in lowered or "skipping" in lowered:
        terminal_warning(text)
    else:
        terminal_print(text)


def _clean_task_description(message):
    text = message.strip().lstrip("-").strip()
    lowered = text.lower()
    stage_descriptions = {
        "optimization parameter collection starting": "Reading optimization files and CPU times",
        "charges collection starting": "Extracting atomic charges and spin densities",
        "fmo collection starting": "Extracting frontier orbitals and molecular moments",
        "nbo parameter collection starting": "Extracting NBO charges and bond orders",
        "nmr parameter collection starting": "Extracting NMR shielding tensors",
        "fukui parameter collection starting": "Calculating Fukui descriptors from charge-state files",
        "mlip parameter collection starting": "Parsing MLIP extxyz descriptor files",
        "steric parameter collection starting": "Calculating steric buried volumes",
    }
    for marker, description in stage_descriptions.items():
        if marker in lowered:
            return description
    if "single point energy collection starting" in lowered:
        return "Updating molecule energies from single-point files"

    replacements = [
        "Parameter Collection starting",
        "Energy Collection starting",
        "Collection starting",
    ]
    for value in replacements:
        text = text.replace(value, "").strip(" -")
    return text or "Processing files"


def _clean_completion_message(message):
    return message.strip().lstrip("-").strip()


def initialize_run_log(output_prefix, version, timestamp, reference, argv=None):
    """Create the single run-level .dat file and write run provenance."""
    logger = Logger(build_log_path(output_prefix), verbose=False)
    logger.write_only(
        f"   MOLDSCRIPT v {version} {timestamp} \n   Citation: {reference}\n"
    )
    command_line = " ".join(["python", "-m", "moldscript", *(argv or [])])
    logger.write_only(f"Command line used in MOLDSCRIPT: {command_line}")
    logger.finalize()


# class for logging
class Logger:
    """Simple file logger used to emit the run-level MOLDSCRIPT.dat audit file."""

    def __init__(self, file_path=None, verbose=True):
        global _RUN_LOG_PATH
        self.verbose = verbose
        self.path = Path(file_path) if file_path else None
        self._handle = None
        self.started_new_file = False
        if self.path:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            key = _log_key(self.path)
            self.started_new_file = key not in _LOG_PATHS_INITIALIZED
            mode = "w" if self.started_new_file else "a"
            self._handle = self.path.open(mode, encoding="utf-8")
            _LOG_PATHS_INITIALIZED.add(key)
            _RUN_LOG_PATH = self.path

    @classmethod
    def silent(cls):
        """Return a logger that swallows all output."""
        return cls(file_path=None, verbose=False)

    def write(self, message):
        """Write a message to the log file and route concise output to the terminal."""
        if self._handle:
            self._handle.write(f"{message}\n")
            self._handle.flush()
        if self.verbose:
            _display_log_message(message)

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


def build_log_path(output_prefix: str, module_code: str = None, suffix: str = "dat") -> Path:
    """Return the filesystem path for the single run audit log."""
    prefix = output_prefix or ""
    filename = f"{prefix}MOLDSCRIPT.{suffix}"
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


def cpu_times_seconds(cpu_times):
    return sum(span.total_seconds() for span in _normalise_cpu_spans(cpu_times))


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


def normalize_workers(workers):
    """Return a conservative worker count from CLI/varfile input."""
    try:
        workers = int(workers)
    except (TypeError, ValueError):
        return 1
    return max(1, workers)


def report_job_progress(logger, current, total, last_step=0):
    """Emit the same 5% progress markers used by the terminal progress parser."""
    if total <= 0:
        return last_step
    percent = int((current / total) * 100)
    step = percent // 5
    if step <= last_step:
        return last_step
    message = f"Progress: {step * 5}% ({current}/{total})"
    if logger:
        logger.write(message)
    else:
        terminal.update_progress(current, total)
    return step


def run_file_jobs(jobs, worker, workers=1, logger=None):
    """
    Run independent file jobs sequentially or in a process pool.

    Results are returned in input order so downstream descriptor tables remain
    deterministic even when jobs complete out of order.
    """
    jobs = list(jobs)
    total = len(jobs)
    if total == 0:
        return []

    worker_count = min(normalize_workers(workers), total)
    results = [None] * total
    last_step = 0

    if worker_count == 1:
        for index, job in enumerate(jobs):
            results[index] = worker(job)
            last_step = report_job_progress(logger, index + 1, total, last_step)
        return results

    with ProcessPoolExecutor(max_workers=worker_count) as executor:
        future_to_index = {
            executor.submit(worker, job): index for index, job in enumerate(jobs)
        }
        completed = 0
        for future in as_completed(future_to_index):
            index = future_to_index[future]
            results[index] = future.result()
            completed += 1
            last_step = report_job_progress(logger, completed, total, last_step)
    return results


def molecule_keys(data_dict):
    """Return only initialized molecule entries, excluding bookkeeping keys."""
    return [
        key
        for key, value in data_dict.items()
        if isinstance(value, dict)
        and {"mol", "atom", "bond"}.issubset(value.keys())
    ]


def filename_match_candidates(fullname):
    stem = Path(str(fullname).replace("\\", "/")).name
    if "." in stem:
        stem = stem.rsplit(".", 1)[0]

    candidates = [stem]
    tempname = stem
    while "_" in tempname:
        tempname = tempname.rsplit("_", 1)[0]
        candidates.append(tempname)
    return candidates


def _preview_values(values, max_items=8):
    values = list(values)
    preview = ", ".join(str(value) for value in values[:max_items])
    if len(values) > max_items:
        preview += f", ... ({len(values)} total)"
    return preview or "none"


def filename_match_error(fullname, data_dict, module_name=None):
    module_label = module_name or "module"
    module_code = module_label.lower()
    if module_code == "fukui":
        suffix_flag = "--suffix_fukui_neutral / --suffix_fukui_reduced / --suffix_fukui_oxidized"
        suffix_var = "suffix_fukui_neutral / suffix_fukui_reduced / suffix_fukui_oxidized"
        example_option = suffix_flag
    else:
        suffix_flag = f"--suffix_{module_code}"
        suffix_var = f"suffix_{module_code}"
        example_option = f"{suffix_flag} {module_code}"
    available = molecule_keys(data_dict)
    candidates = filename_match_candidates(fullname)

    return (
        f"x  Could not match {module_label} file key '{fullname}' to an existing molecule.\n"
        f"   Tried keys: {_preview_values(candidates)}\n"
        f"   Existing molecule keys: {_preview_values(available)}\n"
        f"   This usually means {suffix_flag} is missing or incorrect, or the baseline "
        f"optimization keys were created without --suffix_opt.\n"
        f"   Example: if a {module_label} file is named arbr31_wb97xd_{module_code}.log "
        f"and the matching molecule key should be arbr31_wb97xd, pass the matching "
        f"suffix option (for example {example_option}) or set {suffix_var} in a varfile."
    )


def resolve_data_key(fullname, data_dict, module_name=None, logger=None):
    """Resolve a module filename/key to an initialized molecule key."""
    available = set(molecule_keys(data_dict))
    for candidate in filename_match_candidates(fullname):
        if candidate in available:
            return candidate

    message = filename_match_error(fullname, data_dict, module_name=module_name)
    if logger:
        logger.write(message)
    else:
        emit(message, style="bold red")
    raise SystemExit


def _parse_structure_job(job):
    file_name, source_path = job
    try:
        parsed_data = parse_cc_data(file_name, source_path)
        try:
            mol = xyz2mol.xyz2mol(
                parsed_data.atomnos.tolist(),
                parsed_data.atomcoords[-1].tolist(),
                charge=parsed_data.charge,
            )[0]
            smiles = Chem.MolToSmiles(mol)
            warning = None
        except Exception:
            smiles = ""
            warning = "Encountered an issue with the mol embedding. Skipping smiles string."

        cpu_times = parsed_data.metadata.get("cpu_time") if hasattr(parsed_data, "metadata") else None
        return {
            "file_name": file_name,
            "source_path": source_path,
            "smiles": smiles,
            "atomnos": parsed_data.atomnos,
            "bond_length": parsed_data.bond_data_matrix,
            "scfenergy": parsed_data.scfenergies[-1] * eV_to_hartree,
            "cpu_times": cpu_times,
            "warning": warning,
            "error": None,
        }
    except BaseException as exc:
        return {
            "file_name": file_name,
            "source_path": source_path,
            "error": f"Error parsing {file_name}: {exc}",
        }


def initiate_data_dict(data, logger=None, workers=1):
    """
    Initiates a data dictionary to store all the data from the files.
    Progress and informational messages are emitted to `logger` when provided,
    otherwise they fall back to printing to stdout.
    """
    if logger:
        logger.write(f"Initializing data parsing with SMILES and geometry data")
    else:
        emit(f"Initializing data parsing with SMILES and geometry data")

    total = len(data)
    data_dict = {}
    data_dict["CPU_time"] = []

    if total == 0:
        if logger:
            logger.write("No files to process.")
        else:
            emit("No files to process.", style="yellow")
        return data_dict

    jobs = [(file_name, data[file_name]) for file_name in data.keys()]
    for result in run_file_jobs(jobs, _parse_structure_job, workers=workers, logger=logger):
        file_name = result["file_name"]
        if result.get("error"):
            if logger:
                logger.write(f"x  {result['error']}")
            else:
                emit(f"x  {result['error']}", style="bold red")
            raise SystemExit

        data_dict[file_name] = dict()
        data_dict[file_name]["mol"] = dict()
        data_dict[file_name]["atom"] = dict()
        data_dict[file_name]["bond"] = dict()
        if logger:
            logger.write_only(f"o  Initializing structure data from {os.path.basename(file_name)}")
        if result["warning"]:
            if logger:
                logger.write(result["warning"])
            else:
                emit(result["warning"], style="yellow")
        data_dict[file_name]["mol"]["smiles"] = result["smiles"]
        data_dict[file_name]["atom"]["atomnos"] = result["atomnos"]
        data_dict[file_name]["bond"]["bond_length"] = result["bond_length"]
        data_dict[file_name]["mol"]["scfenergy"] = result["scfenergy"]
        record_cpu_time(data_dict, file_name, result["source_path"], result["cpu_times"])

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
            coords = np.asarray(data.atomcoords[-1], dtype=float)
        except:
            coords = np.asarray(data, dtype=float)
        diff = coords[:, None, :] - coords[None, :, :]
        return np.linalg.norm(diff, axis=-1)
def parse_cc_data(file_name, file):
        try:
            parser = cc.io.ccopen(file)
            cc_data = parser.parse()
            setattr(cc_data, "bond_data_matrix", bond_data_matrix(cc_data))
        except Exception as e:
            # raise an informative SystemExit without printing to stdout here;
            # callers will log the error to the module .dat file as appropriate.
            raise SystemExit(f"Error parsing {file_name}: {e}")
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
    return resolve_data_key(fullname, dd)


