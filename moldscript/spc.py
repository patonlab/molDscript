######################################################.
#        This file stores the spc class               #
######################################################.


import sys, os
import time
import cclib as cc
from moldscript.argument_parser import load_variables
from moldscript.utils import eV_to_hartree
import datetime

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
        if len(self.data.keys()) == 0:
            self.args.log.write(f"\nx  Could not find files to obtain information for single point correction")
            self.args.log.finalize()
            sys.exit()
        else:
            self.file_data = self.get_data()

        if create_dat:
            elapsed_time = round(time.time() - start_time_overall, 2)
            self.args.log.write(f"   --- Single Point Energy Collection complete in {elapsed_time} seconds\n")
            self.args.log.finalize()

    def get_data(self):

        self.args.log.write(f"   --- Single Point Energy Collection starting")

        for file_name in self.data.keys():

            spc_data = self.parse_cc_data(file_name, self.data[file_name])

            filename = self.get_filename(file_name)

            try:
                if list(self.data.keys()).index(file_name) == 0:
                    self.args.log.write(f"   Functional used: {spc_data.metadata['functional']}")
                    self.args.log.write(f"   Basis set used: {spc_data.metadata['basis_set']}")
            except:
                pass
            self.args.log.write(f"o  Parsing SPC Energy Data from {os.path.basename(file_name)}")
            self.data_dict[filename]['mol']['scfenergy'] = (
                spc_data.scfenergies[-1] * eV_to_hartree)
            
            if self.data[file_name] in self.data_dict['CPU_time']:
                pass
            else:
                try: 
                    for time in spc_data.metadata["cpu_time"]:
                        self.data_dict[file_name]["CPU_time"] += time  # add cpu time
                    self.data_dict["CPU_time"].append(self.data[file_name])
                except:
                    self.data_dict[file_name]["CPU_time"] = datetime.timedelta(0)  # initialize cpu time
                    for time in spc_data.metadata["cpu_time"]:
                        self.data_dict[file_name]["CPU_time"] += time  # add cpu time
                    self.data_dict['CPU_time'].append(self.data[file_name])
        return self.data_dict

    def parse_cc_data(self, file_name, file):

        try:
            cc_data = cc.io.ccread(file)

        except:
            self.args.log.write(f"\nx  Could not parse {file_name} to obtain spc energy information")
            cc_data = None
        return cc_data
    
    def get_filename(self, fullname):
        flist = list(self.data_dict.keys())
        tempname = fullname
        for i in range(fullname.count("_")):
            try:
                findex = flist.index(tempname)
                keyname = flist[findex]
                return keyname
            except:
                tempname = tempname.rsplit("_", 1)[0]
                print(tempname)
        print(f"Error processing file {fullname}. Ensure consistent naming as described in the docs.")
        raise SystemExit