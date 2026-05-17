######################################################.
#        This file stores the charges class           #
######################################################.


import sys, os
import time
import datetime
from moldscript.argument_parser import load_variables
from moldscript.utils import (
    format_timedelta,
    get_filename,
    initiate_data_dict,
    progress_iter,
    record_cpu_time,
    safe_parse,
)

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
            self.data_dict = initiate_data_dict(self.data, logger=self.args.log)
        if len(self.data.keys()) == 0:
            self.args.log.write(f"\nx  Could not find files to obtain information for charge data")
            sys.exit()
        else:
            self.file_data = self.get_data()

        if create_dat:
            elapsed_time = round(time.time() - start_time_overall, 2)
            self.args.log.write(f"-- Charges Collection complete in {elapsed_time} seconds")

    def get_data(self):

        self.args.log.write(f"-- Charges Collection starting")
        self.module_cpu_seconds = 0.0
        for idx, file_name in progress_iter(self.data.keys(), self.args.log):
            chg_data = safe_parse(self.data[file_name], logger=self.args.log, context="charge data")
            filename = get_filename(file_name, self.data_dict, logger=self.args.log)

            if idx == 1 and chg_data is not None:
                try:
                    self.args.log.write(f"   Functional used: {chg_data.metadata['functional']}")
                    self.args.log.write(f"   Basis set used: {chg_data.metadata['basis_set']}")
                except (AttributeError, KeyError):
                    pass
            self.args.log.write_only(f"o  Parsing Charge Data from {os.path.basename(file_name)}")
            if len(chg_data.atomcharges.keys()) == 1 and 'mulliken' in chg_data.atomcharges:
                self.data_dict[filename]['atom']['mulliken_charge'] = chg_data.atomcharges['mulliken']
            else:
                for i in chg_data.atomcharges.keys():
                    if 'mulliken' not in i and 'sum' not in i:
                        self.data_dict[filename]['atom'][str(i)+'_charge'] = chg_data.atomcharges[i]

            cpu_times = chg_data.metadata.get("cpu_time") if chg_data and hasattr(chg_data, "metadata") else None
            self.module_cpu_seconds += record_cpu_time(self.data_dict, file_name, self.data[file_name], cpu_times)
        module_cpu_td = datetime.timedelta(seconds=self.module_cpu_seconds)
        if self.module_cpu_seconds:
            self.args.log.write(f"-- Charges CPU time: {format_timedelta(module_cpu_td)}")
        return self.data_dict

