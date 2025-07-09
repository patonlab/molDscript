######################################################.
#        This file stores the spc class               #
######################################################.


import sys, os
import time
import cclib as cc
from collections import defaultdict
from moldscript.argument_parser import load_variables
from moldscript.utils import eV_to_hartree, initiate_data_dict
import datetime

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
        if self.data_dict == {}:
            self.data_dict = initiate_data_dict(self.data)
        if len(self.data.keys()) == 0:
            self.args.log.write(f"\nx  Could not find files to obtain information for charge data")
            self.args.log.finalize()
            sys.exit()
        else:
            self.file_data = self.get_data()

        if create_dat:
            elapsed_time = round(time.time() - start_time_overall, 2)
            self.args.log.write(f"-- Charges Collection complete in {elapsed_time} seconds\n")
            self.args.log.finalize()

    def get_data(self):

        print(f"-- Charges Collection starting")
        for file_name in self.data.keys():
            chg_data = self.parse_cc_data(file_name, self.data[file_name])
            filename = self.get_filename(file_name)

            try:
                if list(self.data.keys()).index(file_name) == 0:
                    self.args.log.write(f"   Functional used: {chg_data.metadata['functional']}")
                    self.args.log.write(f"   Basis set used: {chg_data.metadata['basis_set']}")
            except:
                pass
            self.args.log.write(f"o  Parsing Charge Data from {os.path.basename(file_name)}")
            if len(chg_data.atomcharges.keys()) == 1 and 'mulliken' in chg_data.atomcharges:
                self.data_dict[filename]['atom']['mulliken_charge'] = chg_data.atomcharges['mulliken']
            else:
                for i in chg_data.atomcharges.keys():
                    if 'mulliken' not in i and 'sum' not in i:
                        self.data_dict[filename]['atom'][str(i)+'_charge'] = chg_data.atomcharges[i]
                

            if self.data[file_name] in self.data_dict['CPU_time']:
                pass
            else:
                try: 
                    for time in chg_data.metadata["cpu_time"]:
                        self.data_dict[file_name]["CPU_time"] += time  # add cpu time
                    self.data_dict["CPU_time"].append(self.data[file_name])
                except:
                    self.data_dict[file_name]["CPU_time"] = datetime.timedelta(0)  # initialize cpu time
                    for time in chg_data.metadata["cpu_time"]:
                        self.data_dict[file_name]["CPU_time"] += time  # add cpu time
                    self.data_dict['CPU_time'].append(self.data[file_name])
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