######################################################.
#        This file stores the spc class               #
######################################################.


import sys, os
import time
import cclib as cc
from collections import defaultdict
from moldscript.argument_parser import load_variables
from moldscript.utils import eV_to_hartree
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

        self.args.log.write(f"-- Charges Collection starting")
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
        orca6 = False
        with open(file, 'r') as f:
            for line in f:
                if 'Program Version' in line:
                    version = line.split()[2]
                    if version.startswith('6'):
                        orca6 = True
                        break
        print(f"orca6: {orca6}")
        if not orca6:
            parser = cc.io.ccopen(file)

            try:
                cc_data = parser.parse()
            except:
                self.args.log.write(
                    f"\nx  Could not parse {file_name} to obtain spc energy information")
                cc_data = None
        else:
            # try:
            class CCData:
                def __init__(self):
                    self.atomcharges = {}
                    self.metadata = {}
            cc_data = CCData()
            cc_data.metadata["cpu_time"] = ''
            with open(file, 'r') as f:
                for line in f:
                    if 'TOTAL RUN TIME:' in line:
                        time_parts = line.split(':')[1].strip().split()
                        days = int(time_parts[0])
                        hours = int(time_parts[2])
                        minutes = int(time_parts[4])
                        seconds = int(time_parts[6])
                        milliseconds = int(time_parts[8])
                        total_time = datetime.timedelta(days=days, hours=hours, minutes=minutes, seconds=seconds, milliseconds=milliseconds)
                        cc_data.metadata["cpu_time"] = total_time
                    if 'ATOMIC CHARGES' in line:
                        charge_type = line.split()[0].lower()
                        charges_list = []
                        next(f)  # skip the next line
                        for charge_line in f:
                            if charge_line.strip().startswith('Sum') or charge_line.strip() == '':
                                break
                            parts = charge_line.split()
                            print(parts)
                            if parts != []:
                                charge = float(parts[-1])
                                charges_list.append(charge)
                        cc_data.atomcharges[charge_type] = charges_list
            # except:
            #     self.args.log.write(f"\nx  Could not parse {file_name} to obtain spc energy information")
            #     cc_data = None

        return cc_data
    def get_filename(self, fullname):
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
        print(
            f"Error processing file {fullname}. Ensure consistent naming as described in the docs."
        )
        raise SystemExit