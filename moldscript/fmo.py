######################################################.
#        This file stores the spc class               #
######################################################.


import sys, os
import time
import datetime
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
import numpy as np


def _parse_fmo_job(job):
    file_name, source_path, matched_name = job
    try:
        fmo_data = cc.io.ccread(source_path)
        dipole = np.sqrt(np.sum((fmo_data.moments[0] - fmo_data.moments[1]) ** 2, axis=0))
        homo = fmo_data.moenergies[0][fmo_data.homos[0]]
        lumo = fmo_data.moenergies[0][fmo_data.homos[0] + 1]
        softness = lumo - homo
        chemical_potential = (lumo + homo) / 2
        global_electrophilicity = chemical_potential**2 / (2 * softness)
        mol_values = {
            "dipole": dipole,
            "HOMO": homo,
            "LUMO": lumo,
            "HOMO-LUMO_gap": softness,
            "chemical_potential": chemical_potential,
            "global_electrophilicity": global_electrophilicity,
            "global_nucleophilicity": 1 / global_electrophilicity,
        }
        try:
            quadrupole_moments = fmo_data.moments[2]
            quadrupole_matrix = np.array([
                [quadrupole_moments[0], quadrupole_moments[1], quadrupole_moments[2]],
                [quadrupole_moments[1], quadrupole_moments[3], quadrupole_moments[4]],
                [quadrupole_moments[2], quadrupole_moments[4], quadrupole_moments[5]],
            ])
            mol_values["quadrupole_moment_trace"] = np.trace(quadrupole_matrix)
        except Exception:
            mol_values["quadrupole_moment_trace"] = None

        return {
            "file_name": file_name,
            "source_path": source_path,
            "matched_name": matched_name,
            "mol_values": mol_values,
            "metadata": getattr(fmo_data, "metadata", {}),
            "cpu_times": fmo_data.metadata.get("cpu_time") if hasattr(fmo_data, "metadata") else None,
            "error": None,
        }
    except BaseException as exc:
        return {
            "file_name": file_name,
            "source_path": source_path,
            "matched_name": matched_name,
            "error": f"Could not parse {file_name} to obtain FMO and moment information: {exc}",
        }


class fmo:
    """
    Class containing all the functions for the FMO module related to output files
    """

    def __init__(self, data, data_dict, create_dat=True,  **kwargs):

        start_time_overall = time.time()
        # load default and user-specified variables
        self.args = load_variables(kwargs, "FMO", create_dat=create_dat)
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
            self.args.log.write(f"\nx  Could not find files to obtain information for FMO and moment analysis")
            self.args.log.finalize()
            sys.exit()
        else:
            self.file_data = self.get_data()

        if create_dat:
            elapsed_time = round(time.time() - start_time_overall, 2)
            self.args.log.write(f"-- FMO Collection complete in {elapsed_time} seconds\n")
            self.args.log.finalize()

    def get_data(self):

        self.args.log.write(f"-- FMO Collection starting")
        self.module_cpu_seconds = 0.0

        jobs = []
        for file_name in self.data.keys():
            source_path = self.data[file_name]
            matched_name = self.get_filename(file_name)
            jobs.append((file_name, source_path, matched_name))

        for idx, result in enumerate(
            run_file_jobs(jobs, _parse_fmo_job, workers=self.args.workers, logger=self.args.log)
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

            file_name = result["matched_name"]
            self.args.log.write_only(f"o  Parsing FMO and Moment Data from {os.path.basename(file_name)}")
            for key, value in result["mol_values"].items():
                self.data_dict[file_name]["mol"][key] = value

            self.module_cpu_seconds += cpu_times_seconds(result["cpu_times"])
            record_cpu_time(self.data_dict, file_name, result["source_path"], result["cpu_times"])


        return self.data_dict

    def parse_cc_data(self, file_name, file):

        try:
            cc_data = cc.io.ccread(file)
        except:
            self.args.log.write(
                f"\nx  Could not parse {file_name} to obtain spc energy information"
            )
            cc_data = None
        return cc_data


    def get_filename(self, fullname):
        return resolve_data_key(fullname, self.data_dict, module_name="FMO", logger=self.args.log)


