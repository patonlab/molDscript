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
from moldscript.utils import eV_to_hartree, add_cpu_times
from rdkit.Chem import AllChem
import moldscript.xyz2mol as xyz2mol
from rdkit.Chem import rdFreeSASA

class opt:
    """
    Class containing all the functions for the opt module related to Gaussian output files
    """

    def __init__(self, data, data_dict, create_dat=True, ml = True, **kwargs):

        start_time_overall = time.time()
        # load default and user-specified variables
        self.args = load_variables(kwargs, "OPT", create_dat=create_dat)
        self.data = data
        self.data_dict = data_dict
        self.ml = ml

        if len(self.data.keys()) == 0:
            self.args.log.write(
                f"\nx  Could not find files to obtain optimization information. Exiting program"
            )
            self.args.log.finalize()
            sys.exit()
        else:
            self.file_data = self.get_data()

        if create_dat:
            elapsed_time = round(time.time() - start_time_overall, 2)
            try:
                total_cpu = add_cpu_times(self.file_data)
                self.args.log.write(
                    f"\n   QM optimizations complete in {total_cpu} seconds"
                )
            except:
                pass
            self.args.log.write(f"-- Optimization Parameter Collection complete in {elapsed_time} seconds\n")


    def get_data(self):
        mydict = lambda: defaultdict(mydict)

        self.args.log.write(f"-- Optimization Parameter Collection starting")
        self.data_dict["CPU_time"] = []
        test_file = self.data[list(self.data.keys())[0]]
        xtb = False
        with open(test_file, 'r') as f:
            for _ in range(10):  # Read the first 10 lines
                line = f.readline()
                if 'x T B' in line:
                    xtb = True
                    print('- Identified XTB opt file')
                    break

        for i, file_name in enumerate(self.data.keys()):
            nickname = file_name
            self.data_dict[file_name] = dict()
            self.data_dict[file_name]["mol"] = dict()
            self.data_dict[file_name]["atom"] = dict()
            self.data_dict[file_name]["bond"] = dict()
            self.args.log.write(f"o  Parsing Energy & Thermochemistry Data from {os.path.basename(file_name)}")
            if xtb == False:
                # convert log to smiles
                opt_data = self.parse_cc_data(file_name, self.data[file_name])
                mol = xyz2mol.xyz2mol(opt_data.atomnos.tolist(), opt_data.atomcoords[-1].tolist(), charge=opt_data.charge)[0]
                volume = AllChem.ComputeMolVolume(mol)
                smi = Chem.MolToSmiles(mol)
                try:
                    if i == 0:
                        self.args.log.write(f"   Package used: {opt_data.metadata['package']} {opt_data.metadata['package_version']}")
                        self.args.log.write(f"   Functional used: {opt_data.metadata['functional']}")
                        self.args.log.write(f"   Basis set used: {opt_data.metadata['basis_set']}\n")
                except:
                    pass
                file_name = nickname


                self.data_dict[file_name]["mol"]["scfenergy"] = (opt_data.scfenergies[-1] * eV_to_hartree)
                if self.ml == False:
                    self.data_dict[file_name]["mol"]["opt_enthalpy"] = opt_data.enthalpy
                    self.data_dict[file_name]["mol"]["opt_freeenergy"] = opt_data.freeenergy
                self.data_dict[file_name]["mol"]["smiles"] = smi
                self.data_dict[file_name]["atom"]["atomnos"] = opt_data.atomnos
                self.data_dict[file_name]["bond"]["bond_length"] = opt_data.bond_data_matrix
                self.data_dict[file_name]["mol"]["molecular_volume"] = volume
                try: 
                    for time in opt_data.metadata["cpu_time"]:
                        self.data_dict[file_name]["CPU_time"] += time  # add cpu time
                    self.data_dict["CPU_time"].append(self.data[file_name])
                except:
                    self.data_dict[file_name]["CPU_time"] = datetime.timedelta(0)  # initialize cpu time
                    for time in opt_data.metadata["cpu_time"]:
                        self.data_dict[file_name]["CPU_time"] += time  # add cpu time
                    self.data_dict['CPU_time'].append(self.data[file_name])
            elif xtb == True:
                full_name = self.data[file_name]
                # print(full_name)
                # print(file_name)
                with open(full_name, 'r') as f:
                    atom_coords = []
                    atomnos = []
                    parsing_coords = False
                    for line in f:
                        if '          :  net charge                          ' in line:
                            chg = int(line.split()[-2])
                        if '          | TOTAL ENERGY         ' in line:
                            scf_energy = float(line.split()[3])
                            self.data_dict[file_name]["mol"]["scfenergy"] = scf_energy
                        if 'final structure:' in line:
                            parsing_coords = True
                            for _ in range(3):
                                next(f)  # Skip the next 3 lines
                            continue
                        if parsing_coords:
                            parts = line.split()
                            if len(parts) < 4:
                                parsing_coords = False
                                continue
                            atomnos.append(Chem.GetPeriodicTable().GetAtomicNumber(parts[0]))
                            coords = [float(x) for x in parts[1:4]]
                            atom_coords.append(coords)
                        if '* finished run on' in line:
                            for _ in range(4):
                                cpu_time_line = next(f).strip()
                            days, hours, minutes, seconds = map(float, cpu_time_line.split()[2::2])
                            total_seconds = days * 86400 + hours * 3600 + minutes * 60 + seconds
                            self.data_dict[file_name]["CPU_time"] = datetime.timedelta(seconds=total_seconds)
                            self.data_dict["CPU_time"].append(self.data[file_name])
                bond_data_matrix = self.bond_data_matrix(coords)
                self.data_dict[file_name]["bond"]["bond_length"] = bond_data_matrix
                mol = xyz2mol.xyz2mol(atomnos, atom_coords, charge=chg)[0]
                volume = AllChem.ComputeMolVolume(mol)
                smi = Chem.MolToSmiles(mol)
                self.data_dict[file_name]["atom"]["atomnos"] = np.array(atomnos)
                self.data_dict[file_name]["mol"]["smiles"] = smi
                self.data_dict[file_name]["mol"]["molecular_volume"] = volume


        
        return self.data_dict

    def parse_cc_data(self, file_name, file):

        ### parse data
        parser = cc.io.ccopen(file)
        cc_data = parser.parse()

        setattr(cc_data, "bond_data_matrix", self.bond_data_matrix(cc_data))

        return cc_data

    def bond_data_matrix(self, opt_data):
        try:
            coords = opt_data.atomcoords[-1]
        except:
            coords = opt_data
        bond_data_matrix_list = []
        for atom1 in range(len(coords)):
            row = []
            for atom2 in range(len(coords)):
                p1 = np.array(coords[atom1])
                p2 = np.array(coords[atom2])
                squared_dist = np.sum((p1 - p2) ** 2, axis=0)
                dist = np.sqrt(squared_dist)
                row.append(dist)
            bond_data_matrix_list.append(row)
        return bond_data_matrix_list

    def file_base(self, string):

        try:
            int(string[-1])
        except:
            pass
        else:
            return string
        for i in string[::-1]:
            try:
                int(i)
            except:
                pass
            else:
                lastidx = string.rfind(i) + 1

                break
        startidx = string.rfind("/") + 1

        return string[startidx:lastidx]
