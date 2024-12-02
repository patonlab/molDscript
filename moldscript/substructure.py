######################################################.
#      This file stores the substructure class        #
######################################################.


import sys, os
import time
from rdkit.Chem import Draw
from rdkit import Chem
import cclib as cc
from collections import defaultdict
from moldscript.argument_parser import load_variables
import moldscript.xyz2mol as xyz2mol

class substructure:
    """
    Class containing all the functions for the substructure module related to Gaussian output files
    """

    def __init__(self, data, data_dict, substructure, create_dat=True, **kwargs):

        start_time_overall = time.time()
        # load default and user-specified variables
        self.args = load_variables(kwargs, "SUBSTRUCTURE", create_dat=create_dat)
        self.data = data
        self.substructure = substructure
        self.data_dict = data_dict

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
            self.args.log.write(
                f"\n-- Substructure [{self.substructure}] matching complete in {elapsed_time} seconds\n"
            )
            self.args.log.finalize()

    def get_data(self):


        for file_name in self.data.keys():
            index = self.get_mol(self.data[file_name])
            basename = self.file_base(file_name)
            self.args.log.write(
                f"o  Parsing substructure data from {basename}"
            )
            self.data_dict[file_name]['substructure'] = index

        return self.data_dict

    def get_mol(self, file):

        # convert log to smiles
        parser = cc.io.ccopen(file)
        cc_data = parser.parse()
        mol = xyz2mol.xyz2mol(cc_data.atomnos.tolist(), cc_data.atomcoords[-1].tolist(), charge=cc_data.charge)[0]
        print(mol)

        substructure = Chem.MolFromSmarts(self.substructure)
        Draw.MolToImage(substructure, size=(100, 100))
        indexsall = mol.GetSubstructMatches(substructure)
        indexsall = tuple(x + 1 for x in indexsall[0])
        
        return indexsall
    
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