######################################################.
#        This file stores the FUKUI class            #
######################################################.


import sys, os
import time
import datetime
import cclib as cc
from moldscript.argument_parser import load_variables
import numpy as np
from moldscript.utils import (
    eV_to_hartree,
    parse_cc_data,
    record_cpu_time,
    format_timedelta,
    resolve_data_key,
    initiate_data_dict,
    run_file_jobs,
    cpu_times_seconds,
)
import moldscript.xyz2mol as xyz2mol
from rdkit import Chem


def _fukui_npa_data(file, cc_data):
    start_npop = None
    with open(file, "r") as outfile:
        lines = outfile.readlines()
    list_npop = []
    for i, line in enumerate(lines):
        if line.find(" Summary of Natural Population Analysis:") > -1:
            list_npop.append(i + 6)
    if list_npop:
        start_npop = list_npop[0]
    if start_npop is None:
        return None
    nat_charges = []
    end_npop = start_npop + len(cc_data.atomnos)
    for i in range(start_npop, end_npop):
        nat_charges.append(float(lines[i].split()[2]))
    return nat_charges


def _find_first_match(list_a, list_b):
    for element in list_a:
        if element in list_b:
            return element
    return None


def _parse_fukui_state(file_name, source_path):
    if not source_path:
        return None
    cc_data = cc.io.ccread(source_path)
    try:
        natural = _fukui_npa_data(source_path, cc_data)
        if natural is not None:
            cc_data.atomcharges["natural"] = natural
    except Exception:
        pass
    return {
        "source_path": source_path,
        "energy": cc_data.scfenergies[-1] * eV_to_hartree,
        "atomcharges": dict(cc_data.atomcharges),
        "metadata": getattr(cc_data, "metadata", {}),
        "cpu_times": cc_data.metadata.get("cpu_time") if hasattr(cc_data, "metadata") else None,
    }


def _parse_fukui_job(job):
    raw_file_name, matched_name, file_paths = job
    try:
        neutral_data = _parse_fukui_state(raw_file_name, file_paths.get("neutral"))
        oxidized_data = _parse_fukui_state(raw_file_name, file_paths.get("oxidized"))
        reduced_data = _parse_fukui_state(raw_file_name, file_paths.get("reduced"))

        if neutral_data is None or oxidized_data is None or reduced_data is None:
            return {
                "raw_file_name": raw_file_name,
                "matched_name": matched_name,
                "skip": True,
                "error": None,
            }

        chg = _find_first_match(["natural", "hirshfeld", "mulliken"], list(neutral_data["atomcharges"].keys()))
        if chg is None:
            raise ValueError("No compatible charge set found for Fukui calculation")

        reduced_charges = np.array(reduced_data["atomcharges"][chg])
        neutral_charges = np.array(neutral_data["atomcharges"][chg])
        oxidized_charges = np.array(oxidized_data["atomcharges"][chg])
        fplus = -1 * (reduced_charges - neutral_charges)
        fminus = -1 * (neutral_charges - oxidized_charges)
        rad_fukui = (fplus + fminus) / 2

        return {
            "raw_file_name": raw_file_name,
            "matched_name": matched_name,
            "skip": False,
            "charge_type": chg,
            "metadata": neutral_data["metadata"],
            "mol_values": {
                "vertical_ie": oxidized_data["energy"] - neutral_data["energy"],
                "vertical_ea": reduced_data["energy"] - neutral_data["energy"],
            },
            "atom_values": {
                f"oxidized_{chg}_charges": oxidized_charges,
                f"reduced_{chg}_charges": reduced_charges,
                "fplus": fplus,
                "fminus": fminus,
                "frad": rad_fukui,
            },
            "cpu_records": [
                ("neutral", neutral_data["source_path"], neutral_data["cpu_times"]),
                ("reduced", reduced_data["source_path"], reduced_data["cpu_times"]),
                ("oxidized", oxidized_data["source_path"], oxidized_data["cpu_times"]),
            ],
            "error": None,
        }
    except BaseException as exc:
        return {
            "raw_file_name": raw_file_name,
            "matched_name": matched_name,
            "error": f"Could not parse {raw_file_name} to calculate Fukui descriptors: {exc}",
        }


class fukui:
    """
    Class containing all the functions for the fukui module related to Gaussian output files
    """

    def __init__(self, data, data_dicts, create_dat=True, **kwargs):

        start_time_overall = time.time()
        # load default and user-specified variables
        self.args = load_variables(kwargs, "FUKUI", create_dat=create_dat)
        self.data = data
        self.data_dict = data_dicts
        self.module_cpu_seconds = 0.0
        self.module_cpu_seconds = 0.0
        if self.data_dict == {}:
            self.data_dict = self.fukui_data_dict(self.data)

        if len(self.data.keys()) == 0:
            self.args.log.write(f"x  Could not find files to obtain information for calculating Fukui Coefficients\n")
            sys.exit()
        else:
            self.file_data = self.get_data()

        if create_dat:
            elapsed_time = round(time.time() - start_time_overall, 2)
            self.args.log.write(f"-- Fukui Parameter Collection complete in {elapsed_time} seconds\n")
            self.args.log.finalize()

    def get_data(self):

        first = False
        self.args.log.write(f"-- Fukui Parameter Collection starting")
        jobs = []
        for raw_file_name in list(self.data.keys()):
            file_name = resolve_data_key(raw_file_name, self.data_dict, module_name="FUKUI", logger=self.args.log)
            jobs.append((raw_file_name, file_name, self.data[raw_file_name]))

        for result in run_file_jobs(jobs, _parse_fukui_job, workers=self.args.workers, logger=self.args.log):
            raw_file_name = result["raw_file_name"]
            file_name = result["matched_name"]
            if result.get("error"):
                self.args.log.write(f"x  {result['error']}")
                raise SystemExit

            if result.get("skip"):
                self.args.log.write(f"x  Skipping file {raw_file_name} as one either neutral, oxidized or reduced does not exist!")
                continue

            metadata = result["metadata"]
            if first == False:
                try:
                    self.args.log.write(f"   Package used: {metadata['package']} {metadata['package_version']}")
                    self.args.log.write(f"   Functional used: {metadata['functional']}")
                    self.args.log.write(f"   Basis set used: {metadata['basis_set']}\n")
                except: pass
                self.args.log.write(f"   Charges used for FUKUI: {result['charge_type']}")
                first = True

            self.args.log.write_only(f"o  Parsing Fukui data from {raw_file_name}")
            for key, value in result["mol_values"].items():
                self.data_dict[file_name]['mol'][key] = value
            for key, value in result["atom_values"].items():
                self.data_dict[file_name]['atom'][key] = value

            for label, source, cpu_times in result["cpu_records"]:
                self.module_cpu_seconds += cpu_times_seconds(cpu_times)
                record_cpu_time(self.data_dict, file_name, source, cpu_times)
        module_cpu_td = datetime.timedelta(seconds=self.module_cpu_seconds)
        if self.module_cpu_seconds:
            self.args.log.write(f"-- FUKUI CPU time: {format_timedelta(module_cpu_td)}")
        return self.data_dict

    def parse_cc_data(self, file_name, file):

        try:
            cc_data = cc.io.ccread(file)
        except:
            self.args.log.write(f"\nx  Could not parse {file_name} to obtain information for calculating Fukui Coefficients")
            cc_data = None

        try: cc_data.atomcharges["natural"] = self.npa_data(file, cc_data)
        except: pass

        return cc_data

    def npa_data(self, file, cc_data):
        return _fukui_npa_data(file, cc_data)
    def find_first_match(self, list_a, list_b):
        return _find_first_match(list_a, list_b)
    def fukui_data_dict(self,data):
        """
        Initiates a data dictionary to store all the data from the files.
        """
        neutral_files = {
            file_name: states["neutral"]
            for file_name, states in data.items()
            if "neutral" in states
        }
        return initiate_data_dict(
            neutral_files,
            logger=self.args.log,
            workers=self.args.workers,
        )



