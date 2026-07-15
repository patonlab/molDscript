######################################################.
#        This file stores the spc class               #
######################################################.


import sys, os
import time
import datetime
import cclib as cc
from moldscript.argument_parser import load_variables
from moldscript.utils import eV_to_hartree, initiate_data_dict, record_cpu_time, format_timedelta, resolve_data_key

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
            self.data_dict = initiate_data_dict(self.data, logger=self.args.log)
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

        total = len(self.data)
        last_step = 0
        for idx, file_name in enumerate(self.data.keys(), start=1):
            source_path = self.data[file_name]
            percent = int((idx / total) * 100) if total else 100
            step = percent // 5
            if step > last_step:
                for s in range(last_step + 1, step + 1):
                    self.args.log.write(f"Progress: {s * 5}% ({idx}/{total})")
                last_step = step
            spc_data = self.parse_cc_data(file_name, source_path)

            filename = self.get_filename(file_name)

            try:
                if list(self.data.keys()).index(file_name) == 0:
                    self.args.log.write(f"   Functional used: {spc_data.metadata['functional']}")
                    self.args.log.write(f"   Basis set used: {spc_data.metadata['basis_set']}")
            except:
                pass
            self.args.log.write_only(f"o  Parsing SPC Energy Data from {os.path.basename(file_name)}")
            self.data_dict[filename]['mol']['scfenergy'] = (
                spc_data.scfenergies[-1] * eV_to_hartree)

            cpu_times = spc_data.metadata.get("cpu_time") if spc_data and hasattr(spc_data, "metadata") else None
            self.module_cpu_seconds += record_cpu_time(self.data_dict, filename, source_path, cpu_times)
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

