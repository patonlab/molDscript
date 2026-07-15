######################################################.
#        This file stores the spc class               #
######################################################.


import sys, os
import time
import datetime
import cclib as cc
from moldscript.argument_parser import load_variables
from moldscript.utils import (
    eV_to_hartree,
    initiate_data_dict,
    record_cpu_time,
    format_timedelta,
    resolve_data_key,
    run_file_jobs,
    cpu_times_seconds,
)


def _parse_spc_job(job):
    file_name, source_path, matched_name = job
    try:
        spc_data = cc.io.ccread(source_path)
        return {
            "file_name": file_name,
            "source_path": source_path,
            "matched_name": matched_name,
            "scfenergy": spc_data.scfenergies[-1] * eV_to_hartree,
            "metadata": getattr(spc_data, "metadata", {}),
            "cpu_times": spc_data.metadata.get("cpu_time") if hasattr(spc_data, "metadata") else None,
            "error": None,
        }
    except BaseException as exc:
        return {
            "file_name": file_name,
            "source_path": source_path,
            "matched_name": matched_name,
            "error": f"Could not parse {file_name} to obtain spc energy information: {exc}",
        }

class spc:
    """
    Class containing all the functions for the opt module related to Gaussian output files
    """

    def __init__(self, data, data_dict, create_dat=True,  **kwargs):

        start_time_overall = time.time()
        # load default and user-specified variables
        self.args = load_variables(kwargs, "SPC", create_dat=create_dat)
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
            self.args.log.write(f"\nx  Could not find files to obtain information for single point correction")
            self.args.log.finalize()
            sys.exit()
        else:
            self.file_data = self.get_data()

        if create_dat:
            elapsed_time = round(time.time() - start_time_overall, 2)
            module_cpu_td = datetime.timedelta(seconds=self.module_cpu_seconds)
            self.args.log.write("   --- Single Point CPU time: {}".format(format_timedelta(module_cpu_td)))
            self.args.log.write(f"   --- Single Point Energy Collection complete in {elapsed_time} seconds\n")
            self.args.log.finalize()

    def get_data(self):

        self.args.log.write(f"   --- Single Point Energy Collection starting")
        self.module_cpu_seconds = 0.0

        jobs = []
        for file_name in self.data.keys():
            source_path = self.data[file_name]
            filename = self.get_filename(file_name)
            jobs.append((file_name, source_path, filename))

        for idx, result in enumerate(
            run_file_jobs(jobs, _parse_spc_job, workers=self.args.workers, logger=self.args.log)
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
            self.args.log.write_only(f"o  Parsing SPC Energy Data from {os.path.basename(result['file_name'])}")
            filename = result["matched_name"]
            self.data_dict[filename]['mol']['scfenergy'] = result["scfenergy"]

            self.module_cpu_seconds += cpu_times_seconds(result["cpu_times"])
            record_cpu_time(self.data_dict, filename, result["source_path"], result["cpu_times"])
        return self.data_dict

    def parse_cc_data(self, file_name, file):

        try:
            cc_data = cc.io.ccread(file)

        except:
            self.args.log.write(f"\nx  Could not parse {file_name} to obtain spc energy information")
            cc_data = None
        return cc_data

    def get_filename(self, fullname):
        return resolve_data_key(fullname, self.data_dict, module_name="SPC", logger=self.args.log)

