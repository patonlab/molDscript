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
    Class to create a dataframe of parameters.
    """

    def __init__(self, data_dicts, data_type, substructure='', nbo_suffix='SP_NBO', program='gaussian'):
        self.dd = data_dicts
        self.substructure = substructure
        self.program = program
        if data_type == "molecular":
            mol_df = self.get_mol_df()
            self.mol_df = mol_df 
        if data_type == "bond":
            bond_df = self.get_bond_df()
            self.bond_df = bond_df
        if data_type == "atom":
            atom_df = self.get_atom_df()
            self.atom_df = atom_df


    # create a df of bond properties
    def get_bond_df(self):
        bond_csv = 'bond_level.csv'
        
        calced_list = list(self.dd.keys())
        print('\u25A1  AGGREGATING BOND-LEVEL DESCRIPTORS INTO {}'.format(bond_csv))
        
        dict_df = {k: [] for k in ['species', 'atom1', 'atom2', 'wiberg_bond_order']}
        
            
        if all(x in calced_list for x in ['nbo', 'opt']):
            dict_df = {k: [] for k in ['species', 'atom1', 'atom2', 'atom1_type', 'atom2_type', 'wiberg_bond_order', 'bond_length']}      
            nbo_dict = self.dd['nbo'].file_data
            nbo_data=True
        elif 'opt' in calced_list:
            dict_df = {k: [] for k in ['species', 'atom1', 'atom2', 'atom1_type', 'atom2_type', 'bond_length']}
            nbo_data=False
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
            if nbo_data:
                dict = nbo_dict[file_name]
                bo_matrix = dict['bond_order_matrix']
                length_dict = opt_dict[file_name]
                length_matrix = length_dict['bond_length_matrix']
                atom1 = 0
                atoms = length_dict['opt']['atomnos']
                for row_bo in bo_matrix:
                    for bo, index in zip(row_bo, range(len(row_bo))):
                            if index < atom1:
                                dict_df['wiberg_bond_order'].append(bo)
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
                    if self.program=='gaussian':
                        dict_df = {k: [] for k in ['species', 'atom_index', 'atom_type', 'cm5_charge', 'hirshfeld_charge', 'ox_npa_charge', 'ox_cm5_charge', 'ox_hirshfeld_charge', 'red_npa_charge', 'red_cm5_charge', 'red_hirshfeld_charge', 'fukui_plus', 'fukui_minus', 'fukui_rad']}
                    if self.program=='orca':
                        dict_df = {k: [] for k in ['species', 'atom_index', 'atom_type',  'hirshfeld_charge', 'ox_hirshfeld_charge', 'red_hirshfeld_charge', 'fukui_plus', 'fukui_minus', 'fukui_rad']}    
                    for filename in dict.keys():
                        
                        neut_hirsfeld = list(dict[filename]["neutral"]["atomcharges"]["hirsfeld"])
                        atoms = list(dict[filename]['atomnos'])
                        for atom, num in zip(neut_hirsfeld, atoms):
                            dict_df['species'].append(filename)
                            dict_df['atom_index'].append(neut_hirsfeld.index(atom))
                            element = periodictable.elements[num]
                            dict_df['atom_type'].append(element.symbol)
                        if self.program =='gaussian':
                            neut_cm5 = list(dict[filename]["neutral"]["atomcharges"]["cm5"])
                            for atom in neut_cm5:
                                dict_df['cm5_charge'].append(atom)
                            ox_nat = list(dict[filename]["oxidized"]["atomcharges"]["natural"])
                            for atom in ox_nat:
                                dict_df['ox_npa_charge'].append(atom)
                            ox_cm5 = list(dict[filename]["oxidized"]["atomcharges"]["cm5"])
                            for atom in ox_cm5:
                                dict_df['ox_cm5_charge'].append(atom)
                                
                            red_nat = list(dict[filename]["reduced"]["atomcharges"]["natural"])
                            for atom in red_nat:
                                dict_df['red_npa_charge'].append(atom)
                            red_cm5 = list(dict[filename]["reduced"]["atomcharges"]["cm5"])
                            for atom in red_cm5:
                                dict_df['red_cm5_charge'].append(atom)
                            neut_hirsfeld = list(dict[filename]["neutral"]["atomcharges"]["hirsfeld"])
                        for atom in neut_hirsfeld:
                            dict_df['hirshfeld_charge'].append(atom)
                        ox_hirsfeld = list(dict[filename]["oxidized"]["atomcharges"]["hirsfeld"])
                        for atom in ox_hirsfeld:
                            dict_df['ox_hirshfeld_charge'].append(atom)

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

                        if self.program == 'gaussian':
                            properties = ['species', 'smiles', 'energy', 'enthalpy', 'gibbs_energy', 'dipole', 'XX_quadrupole_moment', 'XY_quadrupole_moment', 'XZ_quadrupole_moment', 'YY_quadrupole_moment', 'YZ_quadrupole_moment', 'ZZ_quadrupole_moment', 'HOMO', 'LUMO', 'HOMO-LUMO_gap']
                        if self.program == 'orca':
                            properties = ['species', 'smiles', 'energy', 'enthalpy', 'gibbs_energy', 'dipole', 'HOMO', 'LUMO', 'HOMO-LUMO_gap']
                        if start == False:
                            dict_df = {k: [] for k in properties}
                            start = True
                        dict_df['species'].append(file_name)
                        dict_df['smiles'].append(final_dict['smiles'])
                        dict_df['energy'].append(final_dict['scfenergy'])
                        dict_df['enthalpy'].append(final_dict['enthalpy'])
                        dict_df['gibbs_energy'].append(final_dict['freeenergy'])
                        dict_df['dipole'].append(final_dict['dipole'])
                        if self.program == 'gaussian':
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
                            dict_df['species'].append(file_name)
                            dict_df['SP_ox_energy'].append(oe)
                            dict_df['SP_red_energy'].append(re)
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
                            dict_df['species'].append(file_name)
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
