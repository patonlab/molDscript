######################################################.
#        This file stores the FUKUI class            #
######################################################.


import sys, os
import time
import datetime
import cclib as cc
from moldscript.argument_parser import load_variables
import numpy as np
from moldscript.utils import eV_to_hartree, parse_cc_data, record_cpu_time, format_timedelta
import moldscript.xyz2mol as xyz2mol
from rdkit import Chem

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
        self.module_cpu_seconds = 0.0
        self.module_cpu_seconds = 0.0
        if self.data_dict == {}:
            self.data_dict = self.fukui_data_dict(self.data)

        if len(self.data.keys()) == 0:
            print(f"x  Could not find files to obtain information for calculating Fukui Coefficients\n")
            sys.exit()
        else:
            self.file_data = self.get_data()

        if create_dat:
            elapsed_time = round(time.time() - start_time_overall, 2)
            print(f"-- Fukui Parameter Collection complete in {elapsed_time} seconds\n")

    def get_data(self):

        first = False
        print(f"-- Fukui Parameter Collection starting")
        for file_name in list(self.data.keys()):
            neutral_data, oxidized_data, reduced_data = None, None, None
            if "neutral" in self.data[file_name].keys():
                neutral_data = self.parse_cc_data(file_name, self.data[file_name]["neutral"])
                if first == False:
                    try:
                        print(f"   Package used: {neutral_data.metadata['package']} {neutral_data.metadata['package_version']}")
                        print(f"   Functional used: {neutral_data.metadata['functional']}")
                        print(f"   Basis set used: {neutral_data.metadata['basis_set']}\n")
                    except: pass
            if "oxidized" in self.data[file_name].keys():
                oxidized_data = self.parse_cc_data(file_name, self.data[file_name]["oxidized"])
            if "reduced" in self.data[file_name].keys():
                reduced_data = self.parse_cc_data(file_name, self.data[file_name]["reduced"])
            if neutral_data != None and oxidized_data != None and reduced_data != None:
                print(f"o  Parsing Fukui data from {file_name}")
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
                datasets = (
                    ("neutral", neutral_data, self.data[file_name].get('neutral')),
                    ("reduced", reduced_data, self.data[file_name].get('reduced')),
                    ("oxidized", oxidized_data, self.data[file_name].get('oxidized')),
                )
                for label, dataset, source in datasets:
                    if not dataset:
                        continue
                    cpu_times = dataset.metadata.get('cpu_time') if hasattr(dataset, 'metadata') else None
                    self.module_cpu_seconds += record_cpu_time(self.data_dict, file_name, source, cpu_times)

            except: print(f'!!Could not obtain CPU time for {file_name}, skipping!!')
        module_cpu_td = datetime.timedelta(seconds=self.module_cpu_seconds)
        if self.module_cpu_seconds:
            print(f"-- FUKUI CPU time: {format_timedelta(module_cpu_td)}")
        return self.data_dict

    def parse_cc_data(self, file_name, file):

        try:
            cc_data = cc.io.ccread(file)
        except:
            print(f"\nx  Could not parse {file_name} to obtain information for calculating Fukui Coefficients")
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
    def fukui_data_dict(self,data):
        """
        Initiates a data dictionary to store all the data from the files.
        """
        print(f"Initializing data parsing with SMILES and geometry data")
        data_dict = {}
        for i, file_name in enumerate(data.keys()):
            data_dict[file_name] = dict()
            data_dict[file_name]["mol"] = dict()
            data_dict[file_name]["atom"] = dict()
            data_dict[file_name]["bond"] = dict()
            parsed_data = parse_cc_data(file_name, data[file_name]['neutral'])
            try:
                mol = xyz2mol.xyz2mol(parsed_data.atomnos.tolist(), parsed_data.atomcoords[-1].tolist(), charge=parsed_data.charge)[0]
                smi = Chem.MolToSmiles(mol)
            except:
                print("Encountered an issue with the mol embedding. Skipping smiles string.")
            data_dict[file_name]["mol"]["smiles"] = smi if 'smi' in locals() else ''
            data_dict[file_name]["atom"]["atomnos"] = parsed_data.atomnos
            data_dict[file_name]["bond"]["bond_length"] = parsed_data.bond_data_matrix
            data_dict[file_name]["mol"]["scfenergy"] = (parsed_data.scfenergies[-1] * eV_to_hartree)
        return data_dict



