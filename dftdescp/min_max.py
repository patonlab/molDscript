######################################################.
#        This file stores the min_max class            #
######################################################.

import pandas as pd
import numpy as np
import math
from dftdescp.utils import find_nth

GAS_CONSTANT = 8.3144621  # J / K / mol
J_TO_AU = 4.184 * 627.509541 * 1000.0  # UNIT CONVERSION

class min_max:
    def __init__(self, cut=0.95, temp = 298.15, spc=False,syllables=2 ):
        self.syllables = syllables
        self.cutoff = 1 - cut
        self.temp = temp
        self.spc = spc
        self.get_boltz_dict()
        self.get_mol_min_max()
        self.get_atom_min_max()
        self.get_bond_min_max()


    def get_boltz_dict(self):
        mol_df = pd.read_csv('molecule_level.csv')
        full_names = mol_df['species']
        codenames = []
        for name in full_names:
            ulineidx = find_nth(name, '_', self.syllables)
            codename = name[:ulineidx]
            codenames.append(codename)
        arrnames = np.array(codenames)
        
        done_list = []
        weight_dict = {}
        if self.spc:
            print('\n   --USING SINGLE POINT CORRECTED ENERGIES FOR MIN-MAX-RANGE CUTOFF--')
            
            mol_df['energy'] = mol_df['spc_energy']
            mol_df.drop(columns=['spc_energy'])
        for name in codenames:
            if not name in done_list:
                
                idxes = np.where(arrnames == name) 
                tempdf = mol_df.iloc[idxes]
                if len(tempdf) == 1:
                    filename = list(tempdf['species'])[0]
                    weight_dict[filename] = True
                else:
                    energies = tempdf['energy'] - tempdf['energy'].min()
                    boltz_sum = 0
                    for e, filename in zip(energies, tempdf['species']):
                        boltz_sum += math.exp(-e * J_TO_AU / GAS_CONSTANT  / self.temp)
                    for e, filename in zip(energies, tempdf['species']):
                        weight = math.exp(-e * J_TO_AU / GAS_CONSTANT / self.temp) / boltz_sum
                        if weight >= self.cutoff:
                            weight_dict[filename] = True
                        else:
                            weight_dict[filename] = False
        self.boltz_dict = weight_dict
    def get_mol_min_max(self):
        mol_min_max_csv = 'min_max_molecule_level.csv'
        print('\u25A1  COMPILING MOLECULE-LEVEL MIN, MAX, AND RANGE VALUES INTO {}'.format(mol_min_max_csv))
       
        mol_df = pd.read_csv('molecule_level.csv')
        full_names = mol_df['species']
        
        
        done_list = []
    
        mol_df = pd.read_csv('molecule_level.csv')
        filtered_df = mol_df

        for species in self.boltz_dict.keys():
            if self.boltz_dict[species] == False:        
                filtered_df = filtered_df.drop(filtered_df[filtered_df['species'] == species].index)
        full_names = filtered_df['species']
        codenames = []
        for name in full_names:
            ulineidx = find_nth(name, '_', self.syllables)
            codename = name[:ulineidx]
            codenames.append(codename)
        arrnames = np.array(codenames)
        
        done_list = []
        columns = list(filtered_df.columns)
        columns.remove('species')
        columns.remove('smiles')
        params = ['species', 'smiles']
        for i in columns:
            params.append(i + '_min')
            params.append(i + '_max')
            params.append(i + '_range')

        min_max_dict = {k: [] for k in params}
        for name in codenames:
            if not name in done_list:
                
                idxes = np.where(arrnames == name) 
                tempdf = filtered_df.iloc[idxes]
                if len(tempdf) == 1:
                    min_max_dict['species'].append(name)
                    min_max_dict['smiles'].append(list(tempdf['smiles'])[0])
                    for i in columns:
                        min_max_dict[i + '_min'].append(list(tempdf[i])[0])
                        min_max_dict[i + '_max'].append(list(tempdf[i])[0])
                        min_max_dict[i + '_range'].append(0)

                else:
                    min_max_dict['species'].append(name)
                    min_max_dict['smiles'].append(list(tempdf['smiles'])[0])
                    for i in columns:
                        min = tempdf[i].min()
                        min_max_dict[i + '_min'].append(min)
                        max = tempdf[i].max()
                        min_max_dict[i + '_max'].append(max)
                        range = max - min
                        min_max_dict[i + '_range'].append(range)
        df = pd.DataFrame(min_max_dict)
        df.drop_duplicates(inplace=True)
        df.sort_values(by='species', inplace=True)
        df.to_csv(mol_min_max_csv, index=False)        



    def get_atom_min_max(self):
        atom_min_max_csv = 'min_max_atom_level.csv'
        try:
            atom_df = pd.read_csv('atom_level.csv')
        except:
            print('\u25A1  SKIPPING ATOM LEVEL min-max-range: no atom_level.csv found')
        else:
            print('\u25A1  COMPILING ATOM-LEVEL MIN, MAX, AND RANGE VALUES INTO {}'.format(atom_min_max_csv))

        filtered_df = atom_df
        for species in self.boltz_dict.keys():
            if self.boltz_dict[species] == False:        
                filtered_df = filtered_df.drop(filtered_df[filtered_df['species'] == species].index)
        done_list = []
        columns = list(filtered_df.columns)
        columns.remove('species')
        columns.remove('atom_index')
        columns.remove('atom_type')
        params = ['species', 'atom_index', 'atom_type']
        for i in columns:
            params.append(i + '_min')
            params.append(i + '_max')
            params.append(i + '_range')
        atoms = atom_df['atom_index'].unique()

        min_max_dict = {k: [] for k in params}
        for atom in atoms:
            done_list = []
            spec_atom = filtered_df.loc[filtered_df['atom_index'] == atom]

            full_names = spec_atom['species']
            codenames = []
            for name in full_names:
                ulineidx = find_nth(name, '_', self.syllables)
                codename = name[:ulineidx]
                codenames.append(codename)
            arrnames = np.array(codenames)
            for name in codenames:
                if not name in done_list:
                    
                    idxes = np.where(arrnames == name) 
                    tempdf = spec_atom.iloc[idxes]
                    if len(tempdf) == 1:
                        min_max_dict['species'].append(name)
                        min_max_dict['atom_index'].append(list(tempdf['atom_index'])[0])
                        min_max_dict['atom_type'].append(list(tempdf['atom_type'])[0])
                        for i in columns:
                            min_max_dict[i + '_min'].append(list(tempdf[i])[0])
                            min_max_dict[i + '_max'].append(list(tempdf[i])[0])
                            min_max_dict[i + '_range'].append(0)

                    else:
                        min_max_dict['species'].append(name)
                        min_max_dict['atom_index'].append(list(tempdf['atom_index'])[0])
                        min_max_dict['atom_type'].append(list(tempdf['atom_type'])[0])
                        for i in columns:
                            min = tempdf[i].min()
                            min_max_dict[i + '_min'].append(min)
                            max = tempdf[i].max()
                            min_max_dict[i + '_max'].append(max)
                            range = max - min
                            min_max_dict[i + '_range'].append(range)
        df = pd.DataFrame(min_max_dict)
        df.drop_duplicates(inplace=True)
        df.sort_values(by=['species', 'atom_index'], inplace=True)
        df.to_csv(atom_min_max_csv, index=False)



    def get_bond_min_max(self):
        bond_min_max_csv = 'min_max_bond_level.csv'
        try:
            bond_df = pd.read_csv('bond_level.csv')
        except:
            print('\u25A1  SKIPPING BOND LEVEL min-max-range: no atom_level.csv found')
        else:
            print('\u25A1  COMPILING BOND-LEVEL MIN, MAX, AND RANGE VALUES INTO {}'.format(bond_min_max_csv))

        
        filtered_df = bond_df
        for species in self.boltz_dict.keys():
            if self.boltz_dict[species] == False:        
                filtered_df = filtered_df.drop(filtered_df[filtered_df['species'] == species].index)
        columns = list(filtered_df.columns)
        columns.remove('species')
        columns.remove('atom1')
        columns.remove('atom1_type')
        columns.remove('atom2')
        columns.remove('atom2_type')
        params = ['species', 'atom1', 'atom1_type','atom2', 'atom2_type']
        for i in columns:
            params.append(i + '_min')
            params.append(i + '_max')
            params.append(i + '_range')
        atom1_list = list(filtered_df['atom1'])
        atom2_list = list(filtered_df['atom2'])
        pair_list = []
        for atom1, atom2 in zip(atom1_list, atom2_list):
            pair = (atom1, atom2)
            pair_list.append(pair)
        pair_arr = np.array(pair_list)
        bonds = np.unique(pair_arr, axis=0)

        min_max_dict = {k: [] for k in params}
        for bond in bonds:
            done_list = []
            atom1 = bond[0]
            atom2 = bond[1]
            spec_bond = bond_df.loc[(bond_df['atom1'] == atom1) & (bond_df['atom2'] == atom2)]
            full_names = spec_bond['species']
            codenames = []
            for name in full_names:
                ulineidx = find_nth(name, '_', self.syllables)
                codename = name[:ulineidx]
                codenames.append(codename)
            arrnames = np.array(codenames)
            for name in codenames:
                if not name in done_list:
                    
                    idxes = np.where(arrnames == name) 
                    tempdf = spec_bond.iloc[idxes]
                    if len(tempdf) == 1:
                        min_max_dict['species'].append(name)
                        min_max_dict['atom1'].append(list(tempdf['atom1'])[0])
                        min_max_dict['atom1_type'].append(list(tempdf['atom1_type'])[0])
                        min_max_dict['atom2'].append(list(tempdf['atom2'])[0])
                        min_max_dict['atom2_type'].append(list(tempdf['atom2_type'])[0])
                        for i in columns:
                            min_max_dict[i + '_min'].append(list(tempdf[i])[0])
                            min_max_dict[i + '_max'].append(list(tempdf[i])[0])
                            min_max_dict[i + '_range'].append(0)

                    else:
                        min_max_dict['species'].append(name)
                        min_max_dict['atom1'].append(list(tempdf['atom1'])[0])
                        min_max_dict['atom1_type'].append(list(tempdf['atom1_type'])[0])
                        min_max_dict['atom2'].append(list(tempdf['atom2'])[0])
                        min_max_dict['atom2_type'].append(list(tempdf['atom2_type'])[0])
                        for i in columns:
                            min = tempdf[i].min()
                            min_max_dict[i + '_min'].append(min)
                            max = tempdf[i].max()
                            min_max_dict[i + '_max'].append(max)
                            range = max - min
                            min_max_dict[i + '_range'].append(range)

                        

        df = pd.DataFrame(min_max_dict)
        df.drop_duplicates(inplace=True)
        df.sort_values(by=['species', 'atom1', 'atom2'], inplace=True)
        df.to_csv(bond_min_max_csv, index=False)