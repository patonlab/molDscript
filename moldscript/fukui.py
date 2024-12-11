######################################################.
#        This file stores the FUKUI class            #
######################################################.


import sys, os
import time
import cclib as cc
from collections import defaultdict
from moldscript.argument_parser import load_variables
import numpy as np
from moldscript.utils import eV_to_hartree
import datetime

class fukui:
    """
    Class containing all the functions for the fukui module related to Gaussian output files
    """

    def __init__(self, data, data_dicts, create_dat=True, **kwargs):

        start_time_overall = time.time()
        # load default and user-specified variables
        self.args = load_variables(kwargs, "FUKUI", create_dat=create_dat)
        self.data = data
        self.data_dict = data_dicts

        if len(self.data.keys()) == 0:
            self.args.log.write(f"x  Could not find files to obtain information for calculating Fukui Coefficients\n")
            self.args.log.finalize()
            sys.exit()
        else:
            self.file_data = self.get_data()

        if create_dat:
            elapsed_time = round(time.time() - start_time_overall, 2)
            self.args.log.write(f"-- Fukui Parameter Collection complete in {elapsed_time} seconds\n")
            self.args.log.finalize()

    def get_data(self):

        first = False
        self.args.log.write(f"-- Fukui Parameter Collection starting")
        for file_name in list(self.data.keys()):
            neutral_data, oxidized_data, reduced_data = None, None, None
            if "neutral" in self.data[file_name].keys():
                neutral_data = self.parse_cc_data(file_name, self.data[file_name]["neutral"])
                if first == False:
                    try:                 
                        self.args.log.write(f"   Package used: {neutral_data.metadata['package']} {neutral_data.metadata['package_version']}")
                        self.args.log.write(f"   Functional used: {neutral_data.metadata['functional']}")
                        self.args.log.write(f"   Basis set used: {neutral_data.metadata['basis_set']}\n")
                    except: pass
            if "oxidized" in self.data[file_name].keys():
                oxidized_data = self.parse_cc_data(file_name, self.data[file_name]["oxidized"])
            if "reduced" in self.data[file_name].keys():
                reduced_data = self.parse_cc_data(file_name, self.data[file_name]["reduced"])

            if neutral_data != None and oxidized_data != None and reduced_data != None:
                self.args.log.write(f"o  Parsing Fukui data from {file_name}")
                neut_e = neutral_data.scfenergies[-1] * eV_to_hartree
                red_e = reduced_data.scfenergies[-1] * eV_to_hartree
                ox_e = oxidized_data.scfenergies[-1] * eV_to_hartree
                self.data_dict[file_name]['mol']['vertical_ie'] = ox_e - neut_e
                self.data_dict[file_name]['mol']['vertical_ea'] = red_e - neut_e
                chg = self.find_first_match(["natural", "hirshfeld", "mulliken"], list(neutral_data.atomcharges.keys()))
                if first == False:
                    self.args.log.write(f'   Charges used for FUKUI: {chg}')
                    first = True
                reduced_charges = np.array(reduced_data.atomcharges[chg])
                neutral_charges = np.array(neutral_data.atomcharges[chg])
                oxidized_charges = np.array(oxidized_data.atomcharges[chg])
                self.data_dict[file_name]['atom'][f'oxidized_{chg}_charges'] = oxidized_charges
                self.data_dict[file_name]['atom'][f'reduced_{chg}_charges'] = reduced_charges
                #multiplied by -1 to change from charge density to electron density to meet fukui definition
                fplus = -1 * (reduced_charges - neutral_charges)
                fminus = -1 * (neutral_charges - oxidized_charges)
                rad_fukui = (fplus + fminus)/2
                self.data_dict[file_name]['atom']['fplus'] = (fplus)
                self.data_dict[file_name]['atom']['fminus'] = (fminus)
                self.data_dict[file_name]['atom']['frad'] = (rad_fukui)
            else:
                self.args.log.write(f"x  Skipping file {file_name} as one either neutral, oxidized or reduced does not exist!")
            try:
                for i in ['neutral', 'reduced', 'oxidized']:
                    if self.data[file_name][i] in self.data_dict['CPU_time']:
                        pass
                    else:
                        try: 
                            if i == 'neutral':
                                for time in neutral_data.metadata["cpu_time"]:
                                    self.data_dict[file_name]["CPU_time"] += time  # add cpu time
                            elif i == 'reduced':
                                for time in reduced_data.metadata["cpu_time"]:
                                    self.data_dict[file_name]["CPU_time"] += time
                            elif i == 'oxidized':
                                for time in oxidized_data.metadata["cpu_time"]:
                                    self.data_dict[file_name]["CPU_time"] += time
                            self.data_dict["CPU_time"].append(self.data[file_name][i])

                        except:
                            self.data_dict[file_name]["CPU_time"] = datetime.timedelta(0)  # initialize cpu time
                            if i == 'neutral':
                                for time in neutral_data.metadata["cpu_time"]:
                                    self.data_dict[file_name]["CPU_time"] += time  # add cpu time
                            elif i == 'reduced':
                                for time in reduced_data.metadata["cpu_time"]:
                                    self.data_dict[file_name]["CPU_time"] += time
                            elif i == 'oxidized':
                                for time in oxidized_data.metadata["cpu_time"]:
                                    self.data_dict[file_name]["CPU_time"] += time
                            self.data_dict['CPU_time'].append(self.data[file_name][i])
            except: print(f'!!Could not obtain CPU time for {file_name}, skipping!!')
        return self.data_dict

    def parse_cc_data(self, file_name, file):

        try:
            cc_data = cc.io.ccread(file)
        except:
            self.args.log.write(f"\nx  Could not parse {file_name} to obtain information for calculating Fukui Coefficients")
            cc_data = None

        try: cc_data.atomcharges["natural"] = self.npa_data(file, cc_data)
        except: pass
        
        return cc_data
        
    def npa_data(self, file, cc_data):
        start_npop = None
        outfile = open(file, "r")
        lines = outfile.readlines()
        list_npop = []
        for i, line in enumerate(lines):
            if line.find(" Summary of Natural Population Analysis:") > -1:
                list_npop.append(i + 6)
        start_npop = list_npop[0]
        if start_npop != None:
            nat_charges = []
            end_npop = start_npop + len(cc_data.atomnos)
            for i in range(start_npop, end_npop):
                nat_charges.append(float(lines[i].split()[2]))
        return nat_charges
    def find_first_match(self, list_a, list_b):
        for element in list_a:
            if element in list_b:
                return element
        return None  # No match found


