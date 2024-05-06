######################################################.
#        This file stores the min_max class            #
######################################################.

import pandas as pd
import numpy as np
import math

GAS_CONSTANT = 8.3144621  # J / K / mol
J_TO_AU = 4.184 * 627.509541 * 1000.0  # UNIT CONVERSION

class min_max:
    def __init__(self, cut=0.95, temp = 298.15 ):
        self.cutoff = 1 - cut
        self.temp = temp
        self.get_boltz_dict()
        self.get_mol_min_max()
        self.get_atom_min_max()


    def get_boltz_dict(self):
        mol_df = pd.read_csv('molecule_level.csv')
        full_names = mol_df['species']
        codenames = []
        for name in full_names:
            ulineidx = name.find('_')
            codename = name[:ulineidx]
            codenames.append(codename)
        arrnames = np.array(codenames)
        
        done_list = []
        weight_dict = {}
        
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
        print('\n\u25A1  COMPILING MOLECULE-LEVEL MIN, MAX, AND RANGE VALUES INTO {}'.format(mol_min_max_csv))
       
        mol_df = pd.read_csv('molecule_level.csv')
        full_names = mol_df['species']
        codenames = []
        for name in full_names:
            ulineidx = name.find('_')
            codename = name[:ulineidx]
            codenames.append(codename)
        arrnames = np.array(codenames)
        
        done_list = []
    
        mol_df = pd.read_csv('molecule_level.csv')
        filtered_df = mol_df

        for species in self.boltz_dict.keys():
            if self.boltz_dict[species] == False:        
                filtered_df = filtered_df.drop(filtered_df[filtered_df['species'] == species].index)
        full_names = filtered_df['species']
        codenames = []
        for name in full_names:
            ulineidx = name.find('_')
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
                ulineidx = name.find('_')
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
        df.to_csv(atom_min_max_csv, index=False)