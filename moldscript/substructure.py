######################################################.
#      This file stores the substructure class        #
######################################################.


import sys, os
import time
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
        except ImportError:
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
        # Parse cclib first so we always have an atom count for the
        # "include all atoms" fallback, even if mol-encoding later fails.
        parser = cc.io.ccopen(file)
        cc_data = parser.parse()
        n_atoms = len(cc_data.atomnos.tolist())

        # Try xyz2mol → RDKit; fall back to OpenBabel if importable; else
        # signal failure with mol=None.
        mol = None
        try:
            mol = xyz2mol.xyz2mol(
                cc_data.atomnos.tolist(),
                cc_data.atomcoords[-1].tolist(),
                charge=cc_data.charge,
            )[0]
        except Exception:
            if self.openbabel:
                try:
                    obConversion = ob.OBConversion()
                    obConversion.SetInAndOutFormats("log", "mol")
                    ob_mol = ob.OBMol()
                    obConversion.ReadFile(ob_mol, file)
                    obConversion.WriteFile(ob_mol, file.split(".")[0] + ".mol")
                    obConversion.CloseOutFile()
                    mol = Chem.MolFromMolFile(file.split(".")[0] + ".mol", removeHs=False)
                except Exception:
                    mol = None

        if mol is None:
            self.args.log.write('!Unable to encode structure for substructure match!')
            self.args.log.write('!Including all atoms in the molecule!')
            if not self.openbabel:
                self.args.log.write('!Consider installing OpenBabel for another encoder option!')
            return tuple(x + 1 for x in range(n_atoms))

        substructure = Chem.MolFromSmarts(self.substructure)
        indexsall = mol.GetSubstructMatches(substructure)
        if indexsall == ():
            self.args.log.write('!No substructure match found!')
            self.args.log.write('!Including all atoms in the molecule!')
            return tuple(x + 1 for x in range(n_atoms))
        return tuple(x + 1 for x in indexsall[0])
    
    def file_base(self, string):
        # Returns the substring up to the last trailing digit (used to strip
        # trailing extension/suffix characters). If the string already ends
        # in a digit, return it unchanged.
        try:
            int(string[-1])
        except (ValueError, IndexError):
            pass
        else:
            return string
        for i in string[::-1]:
            try:
                int(i)
            except ValueError:
                pass
            else:
                lastidx = string.rfind(i) + 1
                break
        startidx = string.rfind('/') + 1
        return string[startidx:lastidx]
