######################################################.
#        This file stores the spc class               #
######################################################.


import sys, os
import time
import datetime
import cclib as cc
from moldscript.argument_parser import load_variables
from moldscript.utils import initiate_data_dict, record_cpu_time, format_timedelta

class charges:
    """
    Class containing all the functions for the charges module related to Gaussian output files
    """

    def __init__(self, data, data_dict, create_dat=True,  **kwargs):

        start_time_overall = time.time()
        # load default and user-specified variables
        self.args = load_variables(kwargs, "SPC", create_dat=create_dat)
        self.data = data
        self.data_dict = data_dict
        self.module_cpu_seconds = 0.0
        if self.data_dict == {}:
            self.data_dict = initiate_data_dict(self.data)
        if len(self.data.keys()) == 0:
            print(f"\nx  Could not find files to obtain information for charge data")
            sys.exit()
        else:
            self.file_data = self.get_data()

        if create_dat:
            elapsed_time = round(time.time() - start_time_overall, 2)
            print(f"-- Charges Collection complete in {elapsed_time} seconds")

    def get_data(self):

        print(f"-- Charges Collection starting")
        self.module_cpu_seconds = 0.0
        for file_name in self.data.keys():
            chg_data = self.parse_cc_data(file_name, self.data[file_name])
            filename = self.get_filename(file_name)

            try:
                if list(self.data.keys()).index(file_name) == 0:
                    print(f"   Functional used: {chg_data.metadata['functional']}")
                    print(f"   Basis set used: {chg_data.metadata['basis_set']}")
            except:
                pass
            print(f"o  Parsing Charge Data from {os.path.basename(file_name)}")
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
            print(f"-- Charges CPU time: {format_timedelta(module_cpu_td)}")
        return self.data_dict

    def parse_cc_data(self, file_name, file):

        parser = cc.io.ccopen(file)

        try:
            cc_data = parser.parse()
        except:
            self.args.log.write(
                f"\nx  Could not parse {file_name} to obtain charge energy information")
            cc_data = None
        return cc_data

    def get_filename(self, fullname):
        try:
            flist = list(self.data_dict.keys())
            tempname = fullname
            for i in range(fullname.count("_")+1):
                try:
                    findex = flist.index(tempname)
                    keyname = flist[findex]
                    return keyname
                except:
                    tempname = tempname.rsplit("_", 1)[0]
                    print(tempname)

        except:
            self.args.log.write('Issue matching one of your filenames, make sure you have a charge file for each opt file')
            raise SystemExit

