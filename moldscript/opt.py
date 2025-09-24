######################################################.
#        This file stores the OPT class               #
######################################################.

import datetime
import sys, os
import time
import cclib as cc
from collections import defaultdict
from moldscript.argument_parser import load_variables
import numpy as np
from rdkit import Chem
from moldscript.utils import eV_to_hartree, initiate_data_dict, parse_cc_data, record_cpu_time, format_timedelta
import moldscript.xyz2mol as xyz2mol

class opt:
    """
    Class containing all the functions for the opt module related to Gaussian output files
    """

    def __init__(self, data, data_dict, create_dat=True, **kwargs):
        start_time_overall = time.time()
        # load default and user-specified variables
        self.args = load_variables(kwargs, "OPT", create_dat=create_dat)
        self.data = data
        self.data_dict = data_dict
        self.module_cpu_seconds = 0.0
        if self.data_dict == {}:
            self.data_dict = initiate_data_dict(self.data)
        if len(self.data.keys()) == 0:
            print(
                f"\nx  Could not find files to obtain optimization information. Exiting program"
            )
            sys.exit()
        else:
            self.file_data = self.get_data()

        if create_dat:
            elapsed_time = round(time.time() - start_time_overall, 2)
            module_cpu_td = datetime.timedelta(seconds=self.module_cpu_seconds)
            if self.module_cpu_seconds:
                self.args.log.write(f"\n   QM optimizations CPU time: {format_timedelta(module_cpu_td)}")
            self.args.log.write(f"-- Optimization Parameter Collection complete in {elapsed_time} seconds\n")


    def get_data(self):
        mydict = lambda: defaultdict(mydict)

        self.args.log.write(f"-- Optimization Parameter Collection starting")
        self.module_cpu_seconds = 0.0
        test_file = self.data[list(self.data.keys())[0]]
        xtb = False
        with open(test_file, 'r') as f:
            for _ in f:  # Read the first 10 lines
                line = f.readline()
                if 'x T B' in line:
                    xtb = True
                    print('- Identified XTB opt file')
                    break

        for i, file_name in enumerate(self.data.keys()):
            self.args.log.write(f"o  Parsing CPU time from {os.path.basename(file_name)}")
            if xtb == False:
                # convert log to smiles
                opt_data = parse_cc_data(file_name, self.data[file_name])
                cpu_times = opt_data.metadata.get('cpu_time') if hasattr(opt_data, 'metadata') else None
                self.module_cpu_seconds += record_cpu_time(self.data_dict, file_name, self.data[file_name], cpu_times)
            elif xtb == True:
                full_name = self.data[file_name]
                with open(full_name, 'r') as f:
                    for line in f:
                        if '* finished run on' in line:
                            for _ in range(4):
                                cpu_time_line = next(f).strip()
                            days, hours, minutes, seconds = map(float, cpu_time_line.split()[2::2])
                            total_seconds = days * 86400 + hours * 3600 + minutes * 60 + seconds
                            cpu_span = datetime.timedelta(seconds=total_seconds)
                            self.module_cpu_seconds += record_cpu_time(self.data_dict, file_name, self.data[file_name], [cpu_span])
        return self.data_dict

