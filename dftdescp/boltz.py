######################################################.
#        This file stores the min_max class            #
######################################################.

import pandas as pd
import numpy as np
import math

GAS_CONSTANT = 8.3144621  # J / K / mol
J_TO_AU = 4.184 * 627.509541 * 1000.0  # UNIT CONVERSION

class boltz:
    def __init__(self, temp = 298.15 ):

        self.temp = temp
        self.mol_boltz()
        self.atom_boltz()
        self.bond_boltz()

    def mol_boltz(self):

        ensemble_mol_csv = 'ensemble_molecule_level.csv'
        print('\u25A1  AVERAGING MOLECULE-LEVEL DESCRIPTORS OVER CONFORMERS INTO {}'.format(ensemble_mol_csv))
       
        mol_df = pd.read_csv('molecule_level.csv')
        full_names = mol_df['species']
        codenames = []
        for name in full_names:
            ulineidx = name.find('_')
            codename = name[:ulineidx]
            codenames.append(codename)
        arrnames = np.array(codenames)
        
        done_list = []
        weighted_df = pd.DataFrame()
        self.weight_dict = {}
        
        for name in codenames:
            if not name in done_list:
                
                idxes = np.where(arrnames == name) 
                tempdf = mol_df.iloc[idxes]
                if len(tempdf) == 1:
                    tempdf['species'] = [name]
                    weighted_df = pd.concat([weighted_df, tempdf])
                    self.weight_dict[name] = 1
                else:
                    energies = tempdf['energy'] - tempdf['energy'].min()
                    columns = list(tempdf.columns)
                    wtrow = {k: [] for k in columns}
                    columns.remove('species')
                    columns.remove('smiles')
                    wtrow['species'] = name
                    smiles = list(tempdf['smiles'])
                    wtrow['smiles'] = smiles[0]
                    boltz_sum = 0.0
                    for e in energies:
                        boltz_sum += math.exp(-e * J_TO_AU / GAS_CONSTANT  / self.temp)
                    weights = []
                    for e in energies:
                        weight = math.exp(-e * J_TO_AU / GAS_CONSTANT / self.temp) / boltz_sum
                        weights.append(weight)
                    self.weight_dict[name] = weights
                    for i in columns:
                        wt_val = 0
                        props = list(tempdf[i])
                        missing_values = []
                        tempweights = weights.copy()
                        for value, idx in zip(props[::-1], range(len(props))[::-1]):
                            if math.isnan(value):
                                missing_values.append(idx)
                        if missing_values != []:

                            for index in missing_values:
                                del props[index]
                                del tempweights[index]
                            tempwt_arr = np.array(tempweights)
                            sum = tempwt_arr.sum()

                            tempweights = [x / sum for x in tempweights]

                        for val, wt in zip(props, tempweights):
                            if math.isnan(val):

                                continue
                            contribution = val * wt
                            wt_val += contribution
                        wtrow[i] = wt_val
                df_row = pd.DataFrame(wtrow, index=[0])
                weighted_df = pd.concat([weighted_df, df_row])
                done_list.append(name)
        boltz_mol_df = weighted_df
        boltz_mol_df.drop_duplicates(inplace=True)
        boltz_mol_df.to_csv(ensemble_mol_csv, index=False)
        mol_df['codenames'] = arrnames
    def atom_boltz(self):
        try:
            atom_df = pd.read_csv('atom_level.csv')
        except:
            print('\u25A1  SKIPPING ATOM LEVEL BOLTZMANN AVERAGING: no atom_level.csv found')
        else:

            ensemble_atom_csv = 'ensemble_atom_level.csv'
            print('\u25A1  AVERAGING ATOM-LEVEL DESCRIPTORS OVER CONFORMERS INTO {}'.format(ensemble_atom_csv))
            atom_df = pd.read_csv('atom_level.csv')
            
            atoms = atom_df['atom_index'].unique()
            weighted_df = pd.DataFrame()
            
            for atom in atoms:
                done_list = []
                spec_atom = atom_df.loc[atom_df['atom_index'] == atom]
                full_names = spec_atom['species']
                codenames = []
                for name in full_names:
                    ulineidx = name.find('_')
                    codename = name[:ulineidx]
                    codenames.append(codename)
                arrnames = np.array(codenames)
                for name in codenames:
                    if not name in done_list:
                        
                        idxes = np.where(arrnames == name) 
                        tempdf = spec_atom.iloc[idxes]
                        if len(tempdf) == 1:
                            tempdf['species'] = [name]

                            weighted_df = pd.concat([weighted_df, tempdf])
                        else:
                            weights = self.weight_dict[name]
                            columns = list(tempdf.columns)
                            wtrow = {k: [] for k in columns}
                            columns.remove('species')
                            columns.remove('atom_type')
                            columns.remove('atom_index')
                            atoms = list(tempdf['atom_type'])
                            wtrow['species'] = name
                            wtrow['atom_type'] = atoms[0]
                            wtrow['atom_index'] = atom
                            for i in columns:
                                wt_val = 0
                                props = list(tempdf[i])
                                missing_values = []
                                tempweights = weights.copy()
                                for value, idx in zip(props[::-1], range(len(props))[::-1]):
                                    if math.isnan(value):
                                        missing_values.append(idx)
                                if missing_values != []:

                                    for index in missing_values:
                                        del props[index]
                                        del tempweights[index]
                                    tempwt_arr = np.array(tempweights)
                                    sum = tempwt_arr.sum()

                                    tempweights = [x / sum for x in tempweights]

                                for val, wt in zip(props, tempweights):
                                    if math.isnan(val):

                                        continue
                                    contribution = val * wt
                                    wt_val += contribution
                                wtrow[i] = wt_val
                        df_row = pd.DataFrame(wtrow, index=[0])
                        weighted_df = pd.concat([weighted_df, df_row])
                    done_list.append(name)

            atom_df = weighted_df
            atom_df.drop_duplicates(inplace=True)
            atom_df.to_csv(ensemble_atom_csv, index=False)
    def bond_boltz(self):
        try:
            bond_df = pd.read_csv('bond_level.csv')
        except:
            print('\u25A1  SKIPPING BOND LEVEL BOLTZMANN AVERAGING: no atom_level.csv found')
        else:

            ensemble_bond_csv = 'ensemble_bond_level.csv'
            print('\u25A1  AVERAGING BOND-LEVEL DESCRIPTORS OVER CONFORMERS INTO {}'.format(ensemble_bond_csv))
            bond_df = pd.read_csv('bond_level.csv')
            atom1_list = list(bond_df['atom1'])
            atom2_list = list(bond_df['atom2'])
            pair_list = []
            for atom1, atom2 in zip(atom1_list, atom2_list):
                pair = (atom1, atom2)
                pair_list.append(pair)
            pair_arr = np.array(pair_list)
            bonds = np.unique(pair_arr, axis=0)
            weighted_df = pd.DataFrame()
            
            for bond in bonds:
                done_list = []
                atom1 = bond[0]
                atom2 = bond[1]
                
                spec_bond = bond_df.loc[(bond_df['atom1'] == atom1) & (bond_df['atom2'] == atom2)]
                full_names = spec_bond['species']
                
                codenames = []
                for name in full_names:
                    ulineidx = name.find('_')
                    codename = name[:ulineidx]
                    codenames.append(codename)
                arrnames = np.array(codenames)
                for name in codenames:
                    if not name in done_list:
                        
                        idxes = np.where(arrnames == name) 
                        
                        tempdf = spec_bond.iloc[idxes]

                        if len(tempdf) == 1:
                            tempdf['species'] = [name]

                            weighted_df = pd.concat([weighted_df, tempdf])
                        else:
                            weights = self.weight_dict[name]
                            columns = list(tempdf.columns)
                            wtrow = {k: [] for k in columns}
                            columns.remove('species')
                            columns.remove('atom1_type')
                            columns.remove('atom1')
                            columns.remove('atom2_type')
                            columns.remove('atom2')
                            atom1_idx = list(tempdf['atom1'])
                            atom1_type = list(tempdf['atom1_type'])

                            atom2_idx = list(tempdf['atom2'])
                            atom2_type = list(tempdf['atom2_type'])
                            wtrow['species'] = name
                            wtrow['atom1_type'] = atom1_type[0]
                            wtrow['atom1'] = atom1_idx[0]
                            wtrow['atom2_type'] = atom2_type[0]
                            wtrow['atom2'] = atom2_idx[0]
                            for i in columns:
                                wt_val = 0
                                props = list(tempdf[i])
                                missing_values = []
                                tempweights = weights.copy()
                                for value, idx in zip(props[::-1], range(len(props))[::-1]):
                                    if math.isnan(value):
                                        missing_values.append(idx)
                                if missing_values != []:

                                    for index in missing_values:
                                        del props[index]
                                        del tempweights[index]
                                    tempwt_arr = np.array(tempweights)
                                    sum = tempwt_arr.sum()

                                    tempweights = [x / sum for x in tempweights]

                                for val, wt in zip(props, tempweights):
                                    if math.isnan(val):

                                        continue
                                    contribution = val * wt
                                    wt_val += contribution
                                wtrow[i] = wt_val
                        df_row = pd.DataFrame(wtrow, index=[0])
                        weighted_df = pd.concat([weighted_df, df_row])
                    done_list.append(name)

            bond_df = weighted_df
            bond_df.drop_duplicates(inplace=True)
            bond_df.to_csv(ensemble_bond_csv, index=False)

