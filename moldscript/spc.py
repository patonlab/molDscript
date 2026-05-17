######################################################.
#        This file stores the spc class               #
######################################################.


import sys, os
import time
import datetime
from moldscript.argument_parser import load_variables
from moldscript.utils import (
    eV_to_hartree,
    format_timedelta,
    get_filename,
    initiate_data_dict,
    progress_iter,
    record_cpu_time,
    safe_parse,
)

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

        for idx, file_name in progress_iter(self.data.keys(), self.args.log):
            spc_data = safe_parse(self.data[file_name], logger=self.args.log, context="spc energy")
            filename = get_filename(file_name, self.data_dict, logger=self.args.log)

            if idx == 1 and spc_data is not None:
                try:
                    self.args.log.write(f"   Functional used: {spc_data.metadata['functional']}")
                    self.args.log.write(f"   Basis set used: {spc_data.metadata['basis_set']}")
                except (AttributeError, KeyError):
                    pass
            self.args.log.write_only(f"o  Parsing SPC Energy Data from {os.path.basename(file_name)}")
            self.data_dict[filename]['mol']['scfenergy'] = (
                spc_data.scfenergies[-1] * eV_to_hartree)

            cpu_times = spc_data.metadata.get("cpu_time") if spc_data and hasattr(spc_data, "metadata") else None
            self.module_cpu_seconds += record_cpu_time(self.data_dict, file_name, self.data[file_name], cpu_times)
        return self.data_dict

