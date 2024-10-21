######################################################.
#        This file stores the FUKUI class            #
######################################################.


import sys, os
import time
import datetime
import cclib as cc
from collections import defaultdict
from moldscript.argument_parser import load_variables
from moldscript.utils import add_cpu_times

class fukui:
    """
    Class containing all the functions for the fukui module related to Gaussian output files
    """

    def __init__(self, data, create_dat=True, **kwargs):

        start_time_overall = time.time()
        # load default and user-specified variables
        self.args = load_variables(kwargs, "FUKUI", create_dat=create_dat)
        self.data = data

        if len(self.data.keys()) == 0:
            self.args.log.write(
                f"x  Could not find files to obtain information for calculating Fukui Coefficients\n"
            )
            self.args.log.finalize()
            sys.exit()
        else:
            self.file_data = self.get_data()

        if create_dat:
            elapsed_time = round(time.time() - start_time_overall, 2)
            try:
                total_cpu = add_cpu_times(self.file_data)
                self.args.log.write(f"\n   Fukui calculations complete in {total_cpu} seconds")
            except: pass
            self.args.log.write(
                f"-- Fukui Parameter Collection complete in {elapsed_time} seconds\n"
            )
            self.args.log.finalize()

    def get_data(self):
        mydict = lambda: defaultdict(mydict)
        file_data = mydict()
        first = False

        for i, file_name in enumerate(self.data.keys()):
            neutral_data, oxidized_data, reduced_data = None, None, None
            if "neutral" in self.data[file_name].keys():
                neutral_data = self.parse_cc_data(
                    file_name, self.data[file_name]["neutral"]
                )
                
            if "oxidized" in self.data[file_name].keys():
                oxidized_data = self.parse_cc_data(
                    file_name, self.data[file_name]["oxidized"]
                )
                
            if "reduced" in self.data[file_name].keys():
                reduced_data = self.parse_cc_data(
                    file_name, self.data[file_name]["reduced"]
                )
            
            if i==0:
                rel_dir = self.data[file_name]["neutral"].split(os.getcwd()+'/')[1].split(file_name)[0]    
                self.args.log.write(
                    f"-- Fukui Parameter Collection from {rel_dir}"
                )
                self.args.log.write(f"   Package used: {neutral_data.metadata['package']} {neutral_data.metadata['package_version']}")
                self.args.log.write(f"   Functional used: {neutral_data.metadata['functional']}")
                self.args.log.write(f"   Basis set used: {neutral_data.metadata['basis_set']}\n")

            if neutral_data != None and oxidized_data != None and reduced_data != None:
                self.args.log.write(
                    f"o  Parsing Fukui data from {file_name}"
                )
                file_data[file_name]["neutral"]["atomcharges"]["natural"] = (
                    neutral_data.atomcharges["natural"]
                )
                file_data[file_name]["neutral"]["atomcharges"]["cm5"] = (
                    neutral_data.atomcharges["cm5"]
                )
                file_data[file_name]["neutral"]["atomcharges"]["hirsfeld"] = (
                    neutral_data.atomcharges["hirsfeld"]
                )
                file_data[file_name]['atomnos'] = neutral_data.atomnos

                file_data[file_name]["oxidized"]["atomcharges"]["natural"] = (
                    oxidized_data.atomcharges["natural"]
                )
                file_data[file_name]["oxidized"]["atomcharges"]["cm5"] = (
                    oxidized_data.atomcharges["cm5"]
                )
                file_data[file_name]["oxidized"]["atomcharges"]["hirsfeld"] = (
                    oxidized_data.atomcharges["hirsfeld"]
                )

                file_data[file_name]["reduced"]["atomcharges"]["natural"] = (
                    reduced_data.atomcharges["natural"]
                )
                file_data[file_name]["reduced"]["atomcharges"]["cm5"] = (
                    reduced_data.atomcharges["cm5"]
                )
                file_data[file_name]["reduced"]["atomcharges"]["hirsfeld"] = (
                    reduced_data.atomcharges["hirsfeld"]
                )
            else:
                self.args.log.write(
                    f"x  Skipping file {file_name} as one either neutral, oxidized or reduced does not exist!"
                )

            file_data[file_name]['cpu_time'] = datetime.timedelta(0) # initialize cpu time
            for time in neutral_data.metadata['cpu_time']:
                file_data[file_name]['cpu_time'] += time # add cpu time from IE
            for time in oxidized_data.metadata['cpu_time']:
                file_data[file_name]['cpu_time'] += time # add cpu time from EA
            for time in reduced_data.metadata['cpu_time']:
                file_data[file_name]['cpu_time'] += time # add cpu time from EA

        return file_data

    def parse_cc_data(self, file_name, file):

        ### parse data
        parser = cc.io.ccopen(file)
        try:
            cc_data = parser.parse()
        except:
            self.args.log.write(
                f"\nx  Could not parse {file_name} to obtain information for calculating Fukui Coefficients"
            )
            cc_data = None

        try: cc_data.atomcharges["natural"] = self.npa_data(file, cc_data)
        except: cc_data.atomcharges["natural"] = None

        if self.args.program=='orca':
            try: cc_data.atomcharges["hirsfeld"], cc_data.atomcharges["cm5"] = self.orca_data(file, cc_data), None
            except: cc_data.atomcharges["hirsfeld"], cc_data.atomcharges["cm5"] = None, None

        if self.args.program=='gaussian':
            try: cc_data.atomcharges["hirsfeld"], cc_data.atomcharges["cm5"] = self.gaussian_data(file, cc_data)
            except: cc_data.atomcharges["hirsfeld"], cc_data.atomcharges["cm5"] = None, None
        
        return cc_data
        
    def npa_data(self, file, cc_data):
        start_npop = None
        outfile = open(file, "r")
        lines = outfile.readlines()
        for i, line in enumerate(lines):
            if line.find(" Summary of Natural Population Analysis:") > -1:
                start_npop = i + 6

        if start_npop != None:
            nat_charges = []
            end_npop = start_npop + len(cc_data.atomnos)
            for i in range(start_npop, end_npop):
                nat_charges.append(float(lines[i].split()[2]))
        return nat_charges
    
    def gaussian_data(self, file, cc_data):
        start_npop = None
        lines = open(file, "r").readlines()
        for i, line in enumerate(lines):
            if (
                line.find(
                    "Hirshfeld charges, spin densities, dipoles, and CM5 charges using IRadAn=      5:"
                )
                > -1
            ):
                start_npop = i + 2
        if start_npop != None:
            hers_charges = []
            cm5_charges = []
            end_npop = start_npop + len(cc_data.atomnos)
            for i in range(start_npop, end_npop):
                hers_charges.append(float(lines[i].split()[2]))
                cm5_charges.append(float(lines[i].split()[-1]))
        return hers_charges, cm5_charges
    
    def orca_data(self, file, cc_data):
        start_npop = None
        lines = open(file, "r").readlines()
        for i, line in enumerate(lines):
            if (
                line.find("HIRSHFELD ANALYSIS")
                > -1
            ):
                start_npop = i + 7
        if start_npop != None:
            hers_charges = []
            end_npop = start_npop + len(cc_data.atomnos)
            for i in range(start_npop, end_npop):
                hers_charges.append(float(lines[i].split()[2]))
        return hers_charges
    

