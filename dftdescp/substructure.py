######################################################.
#      This file stores the substructure class        #
######################################################.


import sys, os
import time

# from openbabel import openbabel as ob
from rdkit import Chem
from collections import defaultdict
from argument_parser import load_variables


class substructure:
    """
    Class containing all the functions for the substructure module related to Gaussian output files
    """

    def __init__(self, data, substructure, create_dat=True, **kwargs):

        start_time_overall = time.time()
        # load default and user-specified variables
        self.args = load_variables(kwargs, "SUBSTRUCTURE", create_dat=create_dat)
        self.data = data
        self.substructure = substructure

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
                f"\nTime Collecting SUBSTRUCTURE data: {elapsed_time} seconds\n"
            )
            self.args.log.finalize()

    def get_data(self):
        mydict = lambda: defaultdict(mydict)
        file_data = mydict()

        for file_name in self.data.keys():
            index = ()
            index = self.get_mol(self.data[file_name])
            self.args.log.write(
                f"Finding substructe information for data from {file_name}\n"
            )
            file_data[file_name][self.substructure]["index"] = index
        return file_data

    def get_mol(self, file):
        # obabel convert
        obConversion = ob.OBConversion()
        obConversion.SetInAndOutFormats("log", "mol")
        ob_mol = ob.OBMol()
        obConversion.ReadFile(ob_mol, file)
        obConversion.WriteFile(ob_mol, file.split(".")[0] + ".mol")
        obConversion.CloseOutFile()
        mol = Chem.MolFromMolFile(file.split(".")[0] + ".mol", removeHs=False)

        substructure = Chem.MolFromSmarts(self.substructure)

        indexsall = mol.GetSubstructMatches(substructure)
        os.remove(file.split(".")[0] + ".mol")

        return indexsall
