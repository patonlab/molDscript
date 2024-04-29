######################################################.
#        This file stores the get_df class            #
######################################################.


import pandas as pd
import numpy as np
import scipy.constants as sc
pd.options.mode.chained_assignment = None
import math
import periodictable

GAS_CONSTANT = 8.3144621  # J / K / mol
J_TO_AU = 4.184 * 627.509541 * 1000.0  # UNIT CONVERSION


class get_df:
    """
    Class to create a dataframe of parameters.
    """

    def __init__(self, data_dicts, data_type, substructure='', nbo_suffix='SP_NBO', temp=298.15):
        self.dd = data_dicts
        self.substructure = substructure
        self.temp = temp
        if data_type == "molecular":
            mol_df = self.get_mol_df()
            self.mol_df = mol_df 
        if data_type == "bond":
            bond_df = self.get_bond_df()
            self.bond_df = bond_df
        if data_type == "atom":
            atom_df = self.get_atom_df()
            self.atom_df = atom_df
        if data_type == "boltzmann":
            self.boltzmann_weight()

    # create a df of bond properties
    def get_bond_df(self):
        bond_csv = 'bond_level.csv'
        calced_list = list(self.dd.keys())
        print('\u25A1  AGGREGATING BOND-LEVEL DESCRIPTORS INTO {}'.format(bond_csv))
        
        dict_df = {k: [] for k in ['species', 'atom1', 'atom2', 'wiberg_bond_order']}
        
            
        if all(x in calced_list for x in ['nbo', 'opt']):
            dict_df = {k: [] for k in ['species', 'atom1', 'atom2', 'atom1_type', 'atom2_type', 'wiberg_bond_order', 'bond_length']}
            nbo_dict = self.dd['nbo'].file_data
            opt_dict = self.dd['opt'].file_data

            for file_name in nbo_dict.keys():
                dict = nbo_dict[file_name]
                bo_matrix = dict['bond_order_matrix']
                length_dict = opt_dict[file_name]
                length_matrix = length_dict['bond_length_matrix']
                atom1 = 0
                atoms = length_dict['opt']['atomnos']

                for row_bo, row_length in zip(bo_matrix, length_matrix):
                    for bo, length, index in zip(row_bo, row_length, range(len(row_bo))):
                        if index < atom1:
                            dict_df['species'].append(file_name)
                            dict_df['atom1'].append(atom1)
                            num = atoms[atom1]
                            element = periodictable.elements[num]
                            dict_df['atom1_type'].append(element.symbol)
                            dict_df['atom2'].append(index)
                            num = atoms[index]
                            element = periodictable.elements[num]
                            dict_df['atom2_type'].append(element.symbol)
                            dict_df['wiberg_bond_order'].append(bo)
                            dict_df['bond_length'].append(length)
                    atom1 +=1
        elif 'opt' in calced_list:
            dict_df = {k: [] for k in ['species', 'atom1', 'atom2',  'atom1_type', 'atom2_type', 'bond_length']}
            opt_dict = self.dd['opt'].file_data

            for file_name in opt_dict.keys():
                length_dict = opt_dict[file_name]
                atoms = length_dict['opt']['atomnos']
                length_matrix = length_dict['bond_length_matrix']
                atom1 = 0
                for  row_length in length_matrix:
                    for length, index in zip(row_length, range(len(row_length))):
                        if index < atom1:
                            dict_df['species'].append(file_name)
                            dict_df['atom1'].append(atom1)
                            num = atoms[atom1]
                            element = periodictable.elements[num]
                            dict_df['atom1_type'].append(element.symbol)
                            dict_df['atom2'].append(index)
                            num = atoms[index]
                            element = periodictable.elements[num]
                            dict_df['atom2_type'].append(element.symbol)
                            dict_df['bond_length'].append(length)
                    atom1 +=1
       
        bond_df = pd.DataFrame(dict_df)
        if self.substructure != '':
            print('   (Filtered by user-defined substructure)')
            final_df = pd.DataFrame()
            for filename in self.substructure.keys():
                dict = self.substructure[filename]
                struc = list(dict.keys())[0]
                filename = self.file_base(filename)
                temp_df = bond_df.loc[bond_df['species'] == filename]
                try:
                    indexes = list(dict[struc]['index'][0])
                except:
                    print(f'{filename}: substructure not found, omitting from final bond df\n')
                    continue
                filtered_df = temp_df.loc[temp_df['atom1'].isin(indexes) | temp_df['atom2'].isin(indexes)]
                final_df = pd.concat([final_df, filtered_df])
                
            final_df.to_csv(bond_csv, index=False)
            return bond_df
        bond_df.to_csv(bond_csv, index=False)
        # print("Saved bond properties to 'bond_df.csv'")
        return bond_df
            

    # create a df of atom properties
    def get_atom_df(self):
        atom_csv = 'atom_level.csv'
        print('\u25A1  AGGREGATING ATOM-LEVEL DESCRIPTORS INTO {}'.format(atom_csv))
        atom_df = pd.DataFrame()
        atom_list = ['nbo', 'nmr', 'fukui']
        calced_list = list(self.dd.keys())
        for category in atom_list:
            if category in calced_list:
                dict = self.dd[category].file_data
                if category == 'nbo':
                    dict_df = {k: [] for k in ['species', 'atom_index', 'atom_type', 'npa_charge', 'wiberg_total']}
                    for filename in dict.keys():
                        
                        charges = list(dict[filename]['charges']['npa'])
                        bond_orders = list(dict[filename]['bond_orders'])
                        #print(len(charges))
                        atoms = list(dict[filename]['atomnos'])
                        for charge, bo, num in zip(charges, bond_orders, atoms):
                            dict_df['wiberg_total'].append(bo)
                            dict_df['npa_charge'].append(charge)
                            dict_df['species'].append(filename)
                            dict_df['atom_index'].append(charges.index(charge))
                            element = periodictable.elements[num]
                            dict_df['atom_type'].append(element.symbol)
                            
                elif category == 'nmr':
                    dict_df = {k: [] for k in ['species', 'atom_index', 'atom_type', 'nmr_shielding']}
                    for filename in dict.keys():
                        
                        shields = list(dict[filename]['nmr_shielding'])
                        atoms = list(dict[filename]['atomnos'])
                        for shield, num in zip(shields, atoms):
                            dict_df['nmr_shielding'].append(shield)
                            dict_df['species'].append(filename)
                            dict_df['atom_index'].append(shields.index(shield))
                            element = periodictable.elements[num]
                            dict_df['atom_type'].append(element.symbol)
                elif category == 'fukui':
                    dict_df = {k: [] for k in ['species', 'atom_index', 'atom_type', 'cm5_charge', 'hirshfeld_charge', 'ox_npa_charge', 'ox_cm5_charge', 'ox_hirshfeld_charge', 'red_npa_charge', 'red_cm5_charge', 'red_hirshfeld_charge', 'fukui_plus', 'fukui_minus', 'fukui_rad']}
                    charges = ['natural', 'cm5', 'hirsfeld']
                    for filename in dict.keys():

                        neut_nat = list(dict[filename]["neutral"]["atomcharges"]["natural"])
                        atoms = list(dict[filename]['atomnos'])
                        for atom, num in zip(neut_nat, atoms):
                            dict_df['species'].append(filename)
                            dict_df['atom_index'].append(neut_nat.index(atom))
                            element = periodictable.elements[num]
                            dict_df['atom_type'].append(element.symbol)
                        neut_cm5 = list(dict[filename]["neutral"]["atomcharges"]["cm5"])
                        for atom in neut_cm5:
                            dict_df['cm5_charge'].append(atom)
                        neut_hirsfeld = list(dict[filename]["neutral"]["atomcharges"]["hirsfeld"])
                        for atom in neut_hirsfeld:
                            dict_df['hirshfeld_charge'].append(atom)
                    
                        ox_nat = list(dict[filename]["oxidized"]["atomcharges"]["natural"])
                        for atom in ox_nat:
                            dict_df['ox_npa_charge'].append(atom)
                        ox_cm5 = list(dict[filename]["oxidized"]["atomcharges"]["cm5"])
                        for atom in ox_cm5:
                            dict_df['ox_cm5_charge'].append(atom)
                        ox_hirsfeld = list(dict[filename]["oxidized"]["atomcharges"]["hirsfeld"])
                        for atom in ox_hirsfeld:
                            dict_df['ox_hirshfeld_charge'].append(atom)

                        red_nat = list(dict[filename]["reduced"]["atomcharges"]["natural"])
                        for atom in red_nat:
                            dict_df['red_npa_charge'].append(atom)
                        red_cm5 = list(dict[filename]["reduced"]["atomcharges"]["cm5"])
                        for atom in red_cm5:
                            dict_df['red_cm5_charge'].append(atom)
                        red_hirsfeld = list(dict[filename]["reduced"]["atomcharges"]["hirsfeld"])
                        for atom in red_hirsfeld:
                            dict_df['red_hirshfeld_charge'].append(atom)
                        for neutral, reduced, oxidized in zip(neut_hirsfeld, red_hirsfeld, ox_hirsfeld):
                            fplus =  reduced- neutral
                            dict_df['fukui_plus'].append(fplus)
                            fminus = neutral- oxidized
                            dict_df['fukui_minus'].append(fminus)
                            frad = fplus - fminus
                            dict_df['fukui_rad'].append(frad)       
                dict_df = pd.DataFrame(dict_df)
                if atom_df.empty:
                    atom_df = dict_df
                else:
                    atom_df = atom_df.merge(dict_df,how='left', on=['species','atom_index', 'atom_type'])
        if self.substructure != '':
            print('   (Filtered by user-defined substructure)')
            final_df = pd.DataFrame()
            for filename in self.substructure.keys():
                dict = self.substructure[filename]
                struc = list(dict.keys())[0]
                basename = self.file_base(filename)
                temp_df = atom_df.loc[atom_df['species'] == basename]
                try:
                    indexes = list(dict[struc]['index'][0])
                except:
                    print(f'{basename}: substructure not found, omitting from final df\n')
                    continue
                filtered_df = temp_df.loc[temp_df['atom'].isin(indexes)]
                final_df = pd.concat([final_df, filtered_df])
                
            final_df.to_csv(atom_csv, index=False)
            return atom_df
        atom_df.to_csv(atom_csv, index=False)
        #print("Saved atom properties to 'atom_df.csv'")
        return atom_df
    
    # create a df of mol properties
    def get_mol_df(self):
        mol_csv = 'molecule_level.csv'
        print('\u25A1  AGGREGATING MOLECULE-LEVEL DESCRIPTORS INTO {}\n'.format(mol_csv))
        mol_df = pd.DataFrame()
        mol_list = ["opt", "sp_ieea", "ad_ieea"]
        calced_list = list(self.dd.keys())
        # go through each of the three options for dictionaries of molecular properties
        for category in mol_list:
            # looks to see if these calcs were done
            if category in calced_list:
                dict = self.dd[category].file_data

                start = False
                if category == 'opt':
                    for file_name in dict.keys():
                        final_dict = dict[file_name]['opt']
                        basename = self.file_base(file_name)

                        properties = ['species', 'smiles', 'energy', 'enthalpy', 'gibbs_energy', 'dipole', 'XX_quadrupole_moment', 'XY_quadrupole_moment', 'XZ_quadrupole_moment', 'YY_quadrupole_moment', 'YZ_quadrupole_moment', 'ZZ_quadrupole_moment', 'HOMO', 'LUMO', 'HOMO-LUMO_gap']
                        if start == False:
                            dict_df = {k: [] for k in properties}
                            start = True
                        dict_df['species'].append(basename)
                        dict_df['smiles'].append(final_dict['smiles'])
                        dict_df['energy'].append(final_dict['scfenergy'])
                        dict_df['enthalpy'].append(final_dict['enthalpy'])
                        dict_df['gibbs_energy'].append(final_dict['freeenergy'])
                        dict_df['dipole'].append(final_dict['dipole'])
                        dict_df['XX_quadrupole_moment'].append(final_dict['XX_quadrupole_moment'])
                        dict_df['XY_quadrupole_moment'].append(final_dict['XY_quadrupole_moment'])  
                        dict_df['XZ_quadrupole_moment'].append(final_dict['XZ_quadrupole_moment'])    
                        dict_df['YY_quadrupole_moment'].append(final_dict['YY_quadrupole_moment']) 
                        dict_df['YZ_quadrupole_moment'].append(final_dict['YZ_quadrupole_moment']) 
                        dict_df['ZZ_quadrupole_moment'].append(final_dict['ZZ_quadrupole_moment']) 
                        dict_df['HOMO'].append(final_dict['HOMO'])
                        dict_df['LUMO'].append(final_dict['LUMO'])
                        dict_df['HOMO-LUMO_gap'].append(final_dict['HOMO-LUMO_gap'])
                elif category == 'sp_ieea':
                     start = False
                     for file_name in dict.keys():
                            basename = self.file_base(file_name)
                            final_dict = dict[file_name]
                            neut_row = mol_df.loc[mol_df['species'] == file_name]
                            neut_e = list(neut_row['energy'])[0]
                            oxidized_e = final_dict['ox']['E']
                            reduced_e = final_dict['red']['E']
                            oe = oxidized_e - neut_e
                            re = reduced_e - neut_e
                            if start == False:
                                dict_df = {k: [] for k in ['species', 'SP_ox_energy','SP_red_energy', 'chemical_hardness', 'global_electrophilicity', 'electronegativity']}
                                start=True
                            dict_df['species'].append(basename)
                            dict_df['SP_ox_energy'].append(oe)
                            dict_df['SP_red_energy'].append(re)
                            ################THIS IS NOT RIGHT BECAUSE NO IE/EA############
                            cp = -1*(oe+re)/2
                            hardness = (oe-re)/2
                            electrophilicity = cp**2 / (2*hardness)
                            electronegativity = -1 * cp
                            dict_df['chemical_hardness'].append(hardness)
                            dict_df['global_electrophilicity'].append(electrophilicity)
                            dict_df['electronegativity'].append(electronegativity)
                elif category == 'ad_ieea':
                     start = False
                     for file_name in dict.keys():
                            basename = self.file_base(file_name)
                            final_dict = dict[file_name]
                            neut_row = mol_df.loc[mol_df['species'] == file_name]
                            neut_e = list(neut_row['energy'])[0]
                            oxidized_e = final_dict['ox']['E']
                            reduced_e = final_dict['red']['E']
                            oe = oxidized_e - neut_e
                            re = reduced_e - neut_e

                            if start == False:
                                dict_df = {k: [] for k in ['species', 'AD_ox_energy','AD_red_energy']}
                                start=True
                            dict_df['species'].append(basename)
                            dict_df['AD_ox_energy'].append(oe)
                            dict_df['AD_red_energy'].append(re)
                dict_df = pd.DataFrame(dict_df)
                if mol_df.empty:
                    mol_df = dict_df
                else:
                    mol_df = mol_df.merge(dict_df,how='left', on='species')
        mol_df.to_csv(mol_csv, index=False)
        # print("Saved molecular properties to 'mol_df.csv'\n\n")
        return mol_df
    
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
    
    def boltzmann_weight(self):
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
        weight_dict = {}
        
        for name in codenames:
            if not name in done_list:
                
                idxes = np.where(arrnames == name) 
                tempdf = mol_df.iloc[idxes]
                if len(tempdf) == 1:
                    tempdf['species'] = [name]
                    weighted_df = pd.concat([weighted_df, tempdf])
                    weight_dict[name] = 1
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
                    weight_dict[name] = weights
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
                            weights = weight_dict[name]
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
                            weights = weight_dict[name]
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













