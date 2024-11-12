######################################################.
#        This file stores the get_df class            #
######################################################.
import pandas as pd
import numpy as np
import scipy.constants as sc
pd.options.mode.chained_assignment = None
import math
import periodictable


class get_df:
    """
    Class to create 3 .csv files containing paramaters
    """

    def __init__(self, data_dicts, substructure="", program="gaussian", volume=False):
        self.dd = data_dicts
        self.substructure = substructure
        mol_df = self.get_mol_df()
        bond_df = self.get_bond_df()
        atom_df = self.get_atom_df()

    def get_mol_df(self):
        mol_csv = "molecule_level.csv"
        print(
            "\n\u25A1  AGGREGATING MOLECULE-LEVEL DESCRIPTORS INTO {}".format(mol_csv)
        )
        filenames = list(self.dd.keys())
        data = self.dd

        for fname in filenames:
            mol_level_data = data[fname]["mol"]
            mol_level_data["filename"] = fname
            tmpdf = pd.DataFrame([mol_level_data])
            try:
                moldf = pd.concat([moldf, tmpdf])
            except:
                moldf = tmpdf
        col = moldf.pop("filename")
        for prop in list(moldf.keys()):
            print(f"\t- {prop}")
        moldf.insert(0, "filename", col)
        moldf.to_csv(mol_csv, index=False)

    def get_bond_df(self):
        bond_csv = "bond_level.csv"
        print("\n\u25A1  AGGREGATING BOND-LEVEL DESCRIPTORS INTO {}".format(bond_csv))
        filenames = list(self.dd.keys())
        data = self.dd
        props = list(data[filenames[0]]["bond"].keys())
        for prop in props:
            print(f"\t- {prop}")
        bonddf = pd.DataFrame()
        for fname in filenames:
            atoms = data[fname]["atom"]["atomnos"]
            filedf = pd.DataFrame()
            for prop in props:
                matrix = np.array(data[fname]["bond"][prop])
                atom1_idx, atom2_idx = np.tril_indices_from(matrix, k=-1)
                values = matrix[atom1_idx, atom2_idx]
                fnames = np.array(fname for i in range(len(values)))
                num1 = atoms[atom1_idx]
                vfunc = np.vectorize(self.get_atom_lab)
                element1 = vfunc(num1)
                num2 = atoms[atom2_idx]
                element2 = vfunc(num2)
                atom1_idx += 1
                atom2_idx += 1
                # Create the DataFrame
                tempdf = pd.DataFrame(
                    {
                        "filename": fnames,
                        "atom1_idx": atom1_idx,
                        "atom1": element1,
                        "atom2_idx": atom2_idx,
                        "atom2": element2,
                        str(prop): values,
                    }
                )
                if not filedf.empty:
                    filedf = pd.merge(
                        filedf,
                        tempdf,
                        on=list(filedf.columns.intersection(tempdf.columns)),
                    )
                else:
                    filedf = tempdf
            if not bonddf.empty:
                bonddf = pd.concat([bonddf, filedf], axis=0)
            else:
                bonddf = filedf
        if self.substructure != '':
            print(f'!!Filtering bond data by user defined substructure: {self.substructure}')
            final_df = pd.DataFrame()
            for filename in filenames:
                idx = self.dd[filename]['substructure']
                temp_df = bonddf.loc[bonddf['filename'] == filename]
                filtered_df = temp_df.loc[temp_df['atom1_idx'].isin(idx) | temp_df['atom2_idx'].isin(idx)]
                final_df = pd.concat([final_df, filtered_df])
            bonddf = final_df    
        bonddf.to_csv(bond_csv, index=False)

    def get_atom_df(self):
        atom_csv = "atom_level.csv"
        print("\n\u25A1  AGGREGATING ATOM-LEVEL DESCRIPTORS INTO {}".format(atom_csv))
        filenames = list(self.dd.keys())
        data = self.dd
        props = list(data[filenames[0]]["atom"].keys())
        props.remove("atomnos")
        for prop in props:
            print(f"\t- {prop}")
        atomdf = pd.DataFrame()
        for fname in filenames:
            atom_level_data = data[fname]["atom"]
            atoms = data[fname]["atom"]["atomnos"]
            tempdic = {}
            for prop in props:
                try:
                    values = atom_level_data[prop]
                except:
                    values = ['' for i in range(len(atoms))]
                tempdic[str(prop)] = values
            fnames = [fname for i in range(len(values))]
            tempdic["filename"] = fnames
            atom_idx = [i + 1 for i in range(len(values))]
            tempdic["atom_index"] = atom_idx
            vfunc = np.vectorize(self.get_atom_lab)
            atypes = np.array(atoms)
            a_labs = vfunc(atypes)
            tempdic["atom_type"] = a_labs
            tempdf = pd.DataFrame(tempdic)
            if not atomdf.empty:
                try:atomdf = pd.concat([atomdf, tempdf])
                except:pass
            else:
                atomdf = tempdf
        col = atomdf.pop("filename")
        atomdf.insert(0, "filename", col)
        if self.substructure != '':
            print(f'!!Filtering atom data by user defined substructure: {self.substructure}\n')
            final_df = pd.DataFrame()
            for filename in filenames:
                idx = self.dd[filename]['substructure']
                temp_df = atomdf.loc[atomdf['filename'] == filename]
                filtered_df = temp_df.loc[temp_df['atom_index'].isin(idx)]
                final_df = pd.concat([final_df, filtered_df])
            atomdf = final_df    
        atomdf.to_csv(atom_csv, index=False)

    def get_atom_lab(self, num):
        label = periodictable.elements[num]
        return label
