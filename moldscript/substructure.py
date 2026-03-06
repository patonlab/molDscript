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
        try:
            from openbabel import openbabel as ob
            self.openbabel = True
        except:
            self.openbabel = False
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
        try:
            parser = cc.io.ccopen(file)
            cc_data = parser.parse()
            mol = xyz2mol.xyz2mol(cc_data.atomnos.tolist(), cc_data.atomcoords[-1].tolist(), charge=cc_data.charge)[0]
        except:
            if self.openbabel:
                try:
                    obConversion = ob.OBConversion()
                    obConversion.SetInAndOutFormats("log", "mol")
                    ob_mol = ob.OBMol()
                    obConversion.ReadFile(ob_mol, file)
                    obConversion.WriteFile(ob_mol, file.split(".")[0] + ".mol")
                    obConversion.CloseOutFile()
                    mol = Chem.MolFromMolFile(file.split(".")[0] + ".mol", removeHs=False)
                except:
                    self.args.log.write('!Unable to encode structure for substructure match!')
                    self.args.log.write('!Including all atoms in the molecule!')
                    whole_mol = tuple(x + 1 for x in range(len(cc_data.atomnos.tolist())))
                    return whole_mol
                else:
                    self.args.log.write('!Unable to encode structure for substructure match!')
                    self.args.log.write('!Including all atoms in the molecule!')
                    self.args.log.write('!Consider installing OpenBabel for another encoder optiion!')
                    whole_mol = tuple(x + 1 for x in range(len(cc_data.atomnos.tolist())))
                    return whole_mol

        substructure = Chem.MolFromSmarts(self.substructure)
        Draw.MolToImage(substructure, size=(100, 100))
        indexsall = mol.GetSubstructMatches(substructure)
        if indexsall == ():
            self.args.log.write('!No substructure match found!')
            self.args.log.write('!Including all atoms in the molecule!')
            indexsall = tuple(x + 1 for x in range(len(cc_data.atomnos.tolist())))
        else:
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
