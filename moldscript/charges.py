######################################################.
#        This file stores the spc class               #
######################################################.


import sys, os
import time
import datetime
import numpy as np
import cclib as cc
from moldscript.argument_parser import load_variables
from moldscript.utils import (
    initiate_data_dict,
    record_cpu_time,
    format_timedelta,
    resolve_data_key,
    run_file_jobs,
    cpu_times_seconds,
)


def _parse_gaussian_mulliken_spins(file):
    spin_tables = []
    with open(file, encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()

    i = 0
    while i < len(lines):
        if lines[i].strip() != "Mulliken charges and spin densities:":
            i += 1
            continue

        table = []
        i += 1
        while i < len(lines):
            parts = lines[i].split()
            if len(parts) == 4 and parts[0].isdigit():
                try:
                    table.append(float(parts[3]))
                except ValueError:
                    break
            elif table:
                break
            i += 1

        if table:
            spin_tables.append(table)
        continue

    if not spin_tables:
        return None
    return np.array(spin_tables[-1])


def _parse_charges_job(job):
    file_name, source_path, matched_name = job
    try:
        parser = cc.io.ccopen(source_path)
        chg_data = parser.parse()
        atom_values = {}
        if len(chg_data.atomcharges.keys()) == 1 and "mulliken" in chg_data.atomcharges:
            atom_values["mulliken_charge"] = chg_data.atomcharges["mulliken"]
        else:
            for charge_type in chg_data.atomcharges.keys():
                if "mulliken" not in charge_type and "sum" not in charge_type:
                    atom_values[f"{charge_type}_charge"] = chg_data.atomcharges[charge_type]

        atom_spins = getattr(chg_data, "atomspins", None)
        if not atom_spins:
            mulliken_spins = _parse_gaussian_mulliken_spins(source_path)
            atom_spins = {"mulliken": mulliken_spins} if mulliken_spins is not None else {}
        for spin_type, spins in atom_spins.items():
            atom_values[f"{spin_type}_spin"] = spins

        return {
            "file_name": file_name,
            "source_path": source_path,
            "matched_name": matched_name,
            "atom_values": atom_values,
            "metadata": getattr(chg_data, "metadata", {}),
            "cpu_times": chg_data.metadata.get("cpu_time") if hasattr(chg_data, "metadata") else None,
            "error": None,
        }
    except BaseException as exc:
        return {
            "file_name": file_name,
            "source_path": source_path,
            "matched_name": matched_name,
            "error": f"Could not parse {file_name} to obtain charge information: {exc}",
        }

class charges:
    """
    Class containing all the functions for the charges module related to Gaussian output files
    """

    def __init__(self, data, data_dict, create_dat=True,  **kwargs):

        start_time_overall = time.time()
        # load default and user-specified variables
        self.args = load_variables(kwargs, "CHARGES", create_dat=create_dat)
        self.data = data
        self.data_dict = data_dict
        self.module_cpu_seconds = 0.0
        if self.data_dict == {}:
            self.data_dict = initiate_data_dict(
                self.data,
                logger=self.args.log,
                workers=self.args.workers,
            )
        if len(self.data.keys()) == 0:
            self.args.log.write(f"\nx  Could not find files to obtain information for charge data")
            sys.exit()
        else:
            self.file_data = self.get_data()

        if create_dat:
            elapsed_time = round(time.time() - start_time_overall, 2)
            self.args.log.write(f"-- Charges Collection complete in {elapsed_time} seconds")
            self.args.log.finalize()

    def get_data(self):

        self.args.log.write(f"-- Charges Collection starting")
        self.module_cpu_seconds = 0.0
        jobs = []
        for file_name in self.data.keys():
            source_path = self.data[file_name]
            filename = self.get_filename(file_name)
            jobs.append((file_name, source_path, filename))

        for idx, result in enumerate(
            run_file_jobs(jobs, _parse_charges_job, workers=self.args.workers, logger=self.args.log)
        ):
            if result.get("error"):
                self.args.log.write(f"\nx  {result['error']}")
                raise SystemExit

            if idx == 0:
                metadata = result["metadata"]
                try:
                    self.args.log.write(f"   Functional used: {metadata['functional']}")
                    self.args.log.write(f"   Basis set used: {metadata['basis_set']}")
                except:
                    pass
            self.args.log.write_only(f"o  Parsing Charge Data from {os.path.basename(result['file_name'])}")
            filename = result["matched_name"]
            for key, value in result["atom_values"].items():
                self.data_dict[filename]['atom'][key] = value

            self.module_cpu_seconds += cpu_times_seconds(result["cpu_times"])
            record_cpu_time(self.data_dict, filename, result["source_path"], result["cpu_times"])
        module_cpu_td = datetime.timedelta(seconds=self.module_cpu_seconds)
        if self.module_cpu_seconds:
            self.args.log.write(f"-- Charges CPU time: {format_timedelta(module_cpu_td)}")
        return self.data_dict

    def parse_cc_data(self, file_name, file):

        parser = cc.io.ccopen(file)

        try:
            cc_data = parser.parse()
        except:
            self.args.log.write_only(
                f"\nx  Could not parse {file_name} to obtain charge energy information")
            cc_data = None
        return cc_data

    def get_atom_spins(self, cc_data, file):
        atom_spins = getattr(cc_data, "atomspins", None) if cc_data is not None else None
        if atom_spins:
            return atom_spins

        mulliken_spins = self.parse_gaussian_mulliken_spins(file)
        if mulliken_spins is None:
            return {}
        return {"mulliken": mulliken_spins}

    @staticmethod
    def parse_gaussian_mulliken_spins(file):
        return _parse_gaussian_mulliken_spins(file)

    def get_filename(self, fullname):
        return resolve_data_key(fullname, self.data_dict, module_name="CHARGES", logger=self.args.log)

