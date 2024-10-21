######################################################.
#        This file stores the OPT class               #
######################################################.


import sys, os
import time
import datetime
import cclib as cc
from collections import defaultdict
from moldscript.argument_parser import load_variables
import numpy as np
from openbabel import openbabel as ob
from rdkit import Chem
from moldscript.utils import eV_to_hartree, add_cpu_times


class opt:
    """
    Class containing all the functions for the opt module related to Gaussian output files
    """

    def __init__(self, data, create_dat=True, **kwargs):

        start_time_overall = time.time()
        # load default and user-specified variables
        self.args = load_variables(kwargs, "OPT", create_dat=create_dat)
        self.data = data

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
                self.args.log.write(f"\n   QM optimizations complete in {total_cpu} seconds")
            except: pass
            self.args.log.write(f"-- Optimization Parameter Collection complete in {elapsed_time} seconds\n")
            self.args.log.finalize()

    def get_data(self):
        mydict = lambda: defaultdict(mydict)
        file_data = mydict()

        self.args.log.write(
                    f"-- Optimization Parameter Collection starting"
                )

        for i, file_name in enumerate(self.data.keys()):
            nickname = file_name

            opt_data = self.parse_cc_data(file_name, self.data[file_name])
            file_name = self.data[file_name]

            obConversion = ob.OBConversion()
            
            if self.args.program == 'gaussian':
                obConversion.SetInAndOutFormats("log", "mol")
                ob_mol = ob.OBMol()
                mol = obConversion.ReadFile(ob_mol, file_name)
                obConversion.WriteFile(ob_mol, file_name.split('.')[0]+'.mol')
                obConversion.CloseOutFile()
                mol = Chem.MolFromMolFile(file_name.split('.')[0]+'.mol', removeHs=False)
                os.remove(file_name.split('.')[0]+'.mol')
                smi = Chem.MolToSmiles(mol)
                
                if i == 0:
                    self.args.log.write(f"   Package used: {opt_data.metadata['package']} {opt_data.metadata['package_version']}")
                    self.args.log.write(f"   Functional used: {opt_data.metadata['functional']}")
                    self.args.log.write(f"   Basis set used: {opt_data.metadata['basis_set']}\n")
            
            if self.args.program == 'orca':
                obConversion.SetInAndOutFormats("out", "mol")
                ob_mol = ob.OBMol()
                mol = obConversion.ReadFile(ob_mol, file_name)
                obConversion.WriteFile(ob_mol, file_name.split('.')[0]+'.mol')
                obConversion.CloseOutFile()
                mol = Chem.MolFromMolFile(file_name.split('.')[0]+'.mol', removeHs=False)
                os.remove(file_name.split('.')[0]+'.mol')
                smi = Chem.MolToSmiles(mol)
            
            file_name = nickname
                
            
            self.args.log.write(f"o  Parsing Energy & Thermochemistry Data from {os.path.basename(file_name)}")
            file_data[file_name]["opt"]["scfenergy"] = (
                opt_data.scfenergies[-1] * eV_to_hartree
            )
            file_data[file_name]["opt"]["enthalpy"] = opt_data.enthalpy
            file_data[file_name]["opt"]["freeenergy"] = opt_data.freeenergy
            file_data[file_name]["opt"]["smiles"] = smi
            file_data[file_name]['opt']['atomnos'] = opt_data.atomnos
        
            file_data[file_name]["bond_length_matrix"] = opt_data.bond_data_matrix
            
            file_data[file_name]["opt"]["dipole"] = np.sqrt(np.sum((opt_data.moments[0]-opt_data.moments[1])**2, axis=0))
            file_data[file_name]["opt"]["HOMO"] = opt_data.moenergies[0][opt_data.homos[0]]

            file_data[file_name]["opt"]["LUMO"] = opt_data.moenergies[0][opt_data.homos[0]+1]
            file_data[file_name]["opt"]["HOMO-LUMO_gap"] = file_data[file_name]["opt"]["LUMO"] - file_data[file_name]["opt"]["HOMO"]
            
            if self.args.program=='gaussian':
                file_data[file_name]["opt"]["XX_quadrupole_moment"] = opt_data.moments[2][0]
                file_data[file_name]["opt"]["XY_quadrupole_moment"] = opt_data.moments[2][1]
                file_data[file_name]["opt"]["XZ_quadrupole_moment"] = opt_data.moments[2][2]
                file_data[file_name]["opt"]["YY_quadrupole_moment"] = opt_data.moments[2][3]
                file_data[file_name]["opt"]["YZ_quadrupole_moment"] = opt_data.moments[2][4]
                file_data[file_name]["opt"]["ZZ_quadrupole_moment"] = opt_data.moments[2][5]

            file_data[file_name]['cpu_time'] = datetime.timedelta(0) # initialize cpu time
            for time in opt_data.metadata['cpu_time']:
                file_data[file_name]['cpu_time'] += time # add cpu time
              
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
        setattr(cc_data, "bond_data_matrix", self.bond_data_matrix(cc_data))

        return cc_data
    
    def bond_data_matrix(self, opt_data):
        coords = opt_data.atomcoords[-1]
        bond_data_matrix_list = []
        for atom1 in range(len(coords)):
            row = []
            for atom2 in range(len(coords)):
                p1 = np.array(coords[atom1])
                p2 = np.array(coords[atom2])
                squared_dist = np.sum((p1-p2)**2, axis=0)
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
        startidx = string.rfind('/') +1

        return string[startidx:lastidx]
