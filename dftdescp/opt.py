######################################################.
#        This file stores the OPT class               #
######################################################.


import sys, os
import time
import cclib as cc
from collections import defaultdict
from dftdescp.argument_parser import load_variables
import numpy as np
import openbabel as ob
from rdkit import Chem

eV_to_hartree = 0.0367493


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
                f"\nx  Could not find files to obtain information optimization"
            )
            self.args.log.finalize()
            sys.exit()
        else:
            self.file_data = self.get_data()

        if create_dat:
            elapsed_time = round(time.time() - start_time_overall, 2)
            self.args.log.write(f"   --- Optimization Parameter Collection complete in {elapsed_time} seconds\n")
            self.args.log.finalize()

    def get_data(self):
        mydict = lambda: defaultdict(mydict)
        file_data = mydict()

        self.args.log.write(
                    f"   --- Optimization Parameter Collection starting"
                )

        for file_name in self.data.keys():
            opt_data = self.parse_cc_data(file_name, self.data[file_name])
            obConversion = ob.OBConversion()
            obConversion.SetInAndOutFormats("log", "mol")
            ob_mol = ob.OBMol()
            mol = obConversion.ReadFile(ob_mol, file_name)
            obConversion.WriteFile(ob_mol, file_name.split('.')[0]+'.mol')
            obConversion.CloseOutFile()
            mol = Chem.MolFromMolFile(file_name.split('.')[0]+'.mol', removeHs=False)
            os.remove(file_name.split('.')[0]+'.mol')
            smi = Chem.MolToSmiles(mol)
            
            if list(self.data.keys()).index(file_name) == 0:
                self.args.log.write(f"   Functional used: {opt_data.metadata['functional']}")
                self.args.log.write(f"   Basis set used: {opt_data.metadata['basis_set']}")
            
            self.args.log.write(f"o  Parsing Energy & Thermochemistry Data from {os.path.basename(file_name)}")
            file_name = self.file_base(file_name)
            file_data[file_name]["opt"]["scfenergy"] = (
                opt_data.scfenergies[-1] * eV_to_hartree
            )
            file_data[file_name]["opt"]["enthalpy"] = opt_data.enthalpy
            file_data[file_name]["opt"]["freeenergy"] = opt_data.freeenergy
            file_data[file_name]["opt"]["smiles"] = smi
            file_data[file_name]['opt']['atomnos'] = opt_data.atomnos



            coords = opt_data.atomcoords[-1]
            bond_data_matrix = []
            for atom1 in range(len(coords)):
                row = []
                for atom2 in range(len(coords)):
                    p1 = np.array(coords[atom1])
                    p2 = np.array(coords[atom2])
                    squared_dist = np.sum((p1-p2)**2, axis=0)
                    dist = np.sqrt(squared_dist)
                    row.append(dist)
                bond_data_matrix.append(row)
            
            file_data[file_name]["bond_length_matrix"] = bond_data_matrix
            moments = opt_data.moments
            com = moments[0]
            xyzdipole = moments[1]
            scalar_dipole = np.sqrt(np.sum((com-xyzdipole)**2, axis=0))
            file_data[file_name]["opt"]["dipole"] = scalar_dipole
            quad_moments = moments[2]
            file_data[file_name]["opt"]["XX_quadrupole_moment"] = quad_moments[0]
            file_data[file_name]["opt"]["XY_quadrupole_moment"] = quad_moments[1]
            file_data[file_name]["opt"]["XZ_quadrupole_moment"] = quad_moments[2]
            file_data[file_name]["opt"]["YY_quadrupole_moment"] = quad_moments[3]
            file_data[file_name]["opt"]["YZ_quadrupole_moment"] = quad_moments[4]
            file_data[file_name]["opt"]["ZZ_quadrupole_moment"] = quad_moments[5]
            homo = opt_data.moenergies[0][opt_data.homos[0]]
            lumo = opt_data.moenergies[0][opt_data.homos[0]+1]
            hl_gap = lumo - homo
            file_data[file_name]["opt"]["HOMO"] = homo
            file_data[file_name]["opt"]["LUMO"] = lumo
            file_data[file_name]["opt"]["HOMO-LUMO_gap"] = hl_gap
                
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

        return cc_data
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
