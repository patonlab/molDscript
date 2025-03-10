######################################################.
#        This file stores the spc class               #
######################################################.


import sys, os
import time
import cclib as cc
from collections import defaultdict
from moldscript.argument_parser import load_variables
from moldscript.utils import eV_to_hartree
import numpy as np
import datetime

class fmo:
    """
    Class containing all the functions for the FMO module related to output files
    """

    def __init__(self, data, data_dict, create_dat=True,  **kwargs):
        
        start_time_overall = time.time()
        # load default and user-specified variables
        self.args = load_variables(kwargs, "FMO", create_dat=create_dat)
        self.data = data
        self.data_dict = data_dict
        if len(self.data.keys()) == 0:
            self.args.log.write(f"\nx  Could not find files to obtain information for FMO and moment analysis")
            self.args.log.finalize()
            sys.exit()
        else:
            self.file_data = self.get_data()

        if create_dat:
            elapsed_time = round(time.time() - start_time_overall, 2)
            self.args.log.write(f"-- FMO Collection complete in {elapsed_time} seconds\n")
            self.args.log.finalize()

    def get_data(self):

        self.args.log.write(f"-- FMO Collection starting")

        for file_name in self.data.keys():
            if self.data[file_name].rsplit('.',1)[1] == 'log':
                self.fmo_program = 'gaussian'
            elif self.data[file_name].rsplit('.', 1)[1] =='out':
                self.fmo_program = 'orca'
            fmo_data = self.parse_cc_data(file_name, self.data[file_name])
            file_name = self.get_filename(file_name)
            try:
                if list(self.data.keys()).index(file_name) == 0:
                    self.args.log.write(f"   Functional used: {fmo_data.metadata['functional']}")
                    self.args.log.write(f"   Basis set used: {fmo_data.metadata['basis_set']}")
            except: pass

            self.args.log.write(f"o  Parsing FMO and Moment Data from {os.path.basename(file_name)}")
            if not self.orca6:
                self.data_dict[file_name]["mol"]["dipole"] = np.sqrt(np.sum((fmo_data.moments[0] - fmo_data.moments[1]) ** 2, axis=0))
                self.data_dict[file_name]["mol"]["HOMO"] = fmo_data.moenergies[0][fmo_data.homos[0]]
                self.data_dict[file_name]["mol"]["LUMO"] = fmo_data.moenergies[0][fmo_data.homos[0] + 1]
            else:
                self.data_dict[file_name]["mol"]["dipole(au)"] = fmo_data.dipole
                self.data_dict[file_name]['mol']['HOMO'] = fmo_data.homo
                self.data_dict[file_name]['mol']['LUMO'] = fmo_data.lumo
  
            softness = self.data_dict[file_name]["mol"]["LUMO"]- self.data_dict[file_name]["mol"]["HOMO"]
            self.data_dict[file_name]["mol"]["HOMO-LUMO_gap"] = (softness)
            chemical_potential = (self.data_dict[file_name]["mol"]["LUMO"] + self.data_dict[file_name]["mol"]["HOMO"]) /2
            self.data_dict[file_name]["mol"]["chemical_potential"] = chemical_potential
            glob_electrophilicity = chemical_potential**2 / (2*softness)
            self.data_dict[file_name]["mol"]["global_electrophilicity"] = glob_electrophilicity
            self.data_dict[file_name]["mol"]["global_nucleophilicity"] = 1/glob_electrophilicity

            try:
                quadrupole_moments = fmo_data.moments[2]
                quadrupole_matrix = np.array([
    [quadrupole_moments[0], quadrupole_moments[1], quadrupole_moments[2]],
    [quadrupole_moments[1], quadrupole_moments[3], quadrupole_moments[4]],
    [quadrupole_moments[2], quadrupole_moments[4], quadrupole_moments[5]]
])
                trace = np.trace(quadrupole_matrix)
                self.data_dict[file_name]["mol"]["quadrupole_moment_trace"] = (trace)
            except: 
                self.data_dict[file_name]["mol"]["quadrupole_moment_trace"] = None
    
            if self.data[file_name] in self.data_dict['CPU_time']:
                pass
            else:
                try: 
                    for time in fmo_data.metadata["cpu_time"]:
                        self.data_dict[file_name]["CPU_time"] += time  # add cpu time
                    self.data_dict["CPU_time"].append(self.data[file_name])
                except:
                    self.data_dict[file_name]["CPU_time"] = datetime.timedelta(0)  # initialize cpu time
                    for time in fmo_data.metadata["cpu_time"]:
                        self.data_dict[file_name]["CPU_time"] += time  # add cpu time
                    self.data_dict['CPU_time'].append(self.data[file_name])


        return self.data_dict

    def parse_cc_data(self, file_name, file):
        self.orca6 = False
        with open(file, 'r') as f:
            for line in f:
                if 'Program Version' in line:
                    version = line.split()[2]
                    if version.startswith('6'):
                        self.orca6 = True
                        break   
        if not self.orca6:

            try:
                cc_data = cc.io.ccread(file)
            except:
                self.args.log.write(
                    f"\nx  Could not parse {file_name} to obtain spc energy information"
                )
                cc_data = None
        else:
            class CCData:
                def __init__(self):
                    self.metadata = {}
                    self.homo = float()
                    self.lumo = float()
                    self.dipole = float()
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
                    if 'ORBITAL ENERGIES' in line:
                        next(f)
                        next(f)
                        next(f)
                        line = next(f)
                        moenergies = []
                        mooccupancies = []
                        proceed = True
                        while proceed==True:
                            if line.startswith('*') or line.strip() == '':
                                proceed = False
                                break
                            parts = line.split()
                            mooccupancies.append(float(parts[1]))
                            moenergies.append(float(parts[3]))
                            line=next(f)
                        for i in range(len(mooccupancies) - 1, -1, -1):
                            if mooccupancies[i] > 0.0:
                                cc_data.homo = moenergies[i]
                                cc_data.lumo = moenergies[i + 1] if i + 1 < len(moenergies) else None
                                break
                    if 'Magnitude (a.u.)' in line:
                        cc_data.dipole = float(line.split(':')[1].strip())
                    
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
        self.args.log.write('Issue matching one of your filenames')
        raise SystemExit