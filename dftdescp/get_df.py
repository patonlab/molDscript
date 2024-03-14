######################################################.
#        This file stores the get_df class            #
######################################################.


import pandas as pd



class get_df:
    """
    Class to create a dataframe of parameters.
    """

    def __init__(self, data_dicts, data_type, substructure=''):
        self.dd = data_dicts
        self.substructure = substructure
        #print(self.dd)
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
        print('\u25A1  AGGREGATING BOND-LEVEL DESCRIPTORS INTO {}'.format(bond_csv))
        nbo_dict = self.dd.file_data
        dict_df = {k: [] for k in ['File', 'atom1', 'atom2', 'bonding']}
        for file_name in nbo_dict.keys():
            dict = nbo_dict[file_name]
            bo_matrix = dict['bond_order_matrix']
            atom1 = 0
            for row in bo_matrix:
                for combination, index in zip(row, range(len(row))):
                    if index < atom1:
                        dict_df['File'].append(file_name)
                        dict_df['atom1'].append(atom1)
                        dict_df['atom2'].append(index)
                        dict_df['bonding'].append(combination)
                atom1 +=1
        bond_df = pd.DataFrame(dict_df)
        if self.substructure != '':
            print('Filtering atomic property df by substructure\n\n')
            final_df = pd.DataFrame()
            for filename in self.substructure.keys():
                dict = self.substructure[filename]
                struc = list(dict.keys())[0]
                filename = self.file_base(filename)
                temp_df = bond_df.loc[bond_df['File'] == filename]
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
                    dict_df = {k: [] for k in ['File', 'atom', 'nbo_charge', 'bond_order']}
                    for filename in dict.keys():
                        
                        charges = list(dict[filename]['charges']['npa'])
                        bond_orders = list(dict[filename]['bond_orders'])
                        #print(len(charges))
                        for charge, bo in zip(charges, bond_orders):
                            dict_df['bond_order'].append(bo)
                            dict_df['nbo_charge'].append(charge)
                            dict_df['File'].append(filename)
                            dict_df['atom'].append(charges.index(charge))
                elif category == 'nmr':
                    dict_df = {k: [] for k in ['File', 'atom', 'nmr_shielding']}
                    for filename in dict.keys():
                        
                        shields = list(dict[filename]['nmr_shielding'])

                        for shield in shields:
                            dict_df['nmr_shielding'].append(shield)
                            dict_df['File'].append(filename)
                            dict_df['atom'].append(shields.index(shield))
                elif category == 'fukui':
                    dict_df = {k: [] for k in ['File', 'atom', 'neut_nat', 'neut_cm5', 'neut_mulliken', 'neut_hirsfeld', 'ox_nat', 'ox_cm5', 'ox_mulliken', 'ox_hirsfeld', 'red_nat', 'red_cm5', 'red_mulliken', 'red_hirsfeld']}
                    charges = ['natural', 'cm5', 'mulliken', 'hirsfeld']
                    for filename in dict.keys():

                        neut_nat = list(dict[filename]["neutral"]["atomcharges"]["natural"])
                        for atom in neut_nat:
                            dict_df['File'].append(filename)
                            dict_df['atom'].append(neut_nat.index(atom))
                            dict_df['neut_nat'].append(atom)
                        neut_cm5 = list(dict[filename]["neutral"]["atomcharges"]["cm5"])
                        for atom in neut_cm5:
                            dict_df['neut_cm5'].append(atom)
                        neut_mulliken = list(dict[filename]["neutral"]["atomcharges"]["mulliken"])
                        for atom in neut_mulliken:
                            dict_df['neut_mulliken'].append(atom)
                        neut_hirsfeld = list(dict[filename]["neutral"]["atomcharges"]["hirsfeld"])
                        for atom in neut_hirsfeld:
                            dict_df['neut_hirsfeld'].append(atom)

                        ox_nat = list(dict[filename]["oxidized"]["atomcharges"]["natural"])
                        for atom in ox_nat:
                            dict_df['ox_nat'].append(atom)
                        ox_cm5 = list(dict[filename]["oxidized"]["atomcharges"]["cm5"])
                        for atom in ox_cm5:
                            dict_df['ox_cm5'].append(atom)
                        ox_mulliken = list(dict[filename]["oxidized"]["atomcharges"]["mulliken"])
                        for atom in ox_mulliken:
                            dict_df['ox_mulliken'].append(atom)
                        ox_hirsfeld = list(dict[filename]["oxidized"]["atomcharges"]["hirsfeld"])
                        for atom in ox_hirsfeld:
                            dict_df['ox_hirsfeld'].append(atom)

                        red_nat = list(dict[filename]["reduced"]["atomcharges"]["natural"])
                        for atom in red_nat:
                            dict_df['red_nat'].append(atom)
                        red_cm5 = list(dict[filename]["reduced"]["atomcharges"]["cm5"])
                        for atom in red_cm5:
                            dict_df['red_cm5'].append(atom)
                        red_mulliken = list(dict[filename]["reduced"]["atomcharges"]["mulliken"])
                        for atom in red_mulliken:
                            dict_df['red_mulliken'].append(atom)
                        red_hirsfeld = list(dict[filename]["reduced"]["atomcharges"]["hirsfeld"])
                        for atom in red_hirsfeld:
                            dict_df['red_hirsfeld'].append(atom)

                dict_df = pd.DataFrame(dict_df)
                if atom_df.empty:
                    atom_df = dict_df
                else:
                    atom_df = atom_df.merge(dict_df,how='left', on=['File','atom'])
        if self.substructure != '':
            print('Filtering atomic property df by substructure\n\n')
            final_df = pd.DataFrame()
            for filename in self.substructure.keys():
                dict = self.substructure[filename]
                struc = list(dict.keys())[0]
                basename = self.file_base(filename)
                temp_df = atom_df.loc[atom_df['File'] == basename]
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
                        for title in dict[file_name].keys():

                            basename = self.file_base(file_name)
                            final_dict = dict[file_name][title]

                            properties = list(final_dict.keys())

                            properties.insert(0, 'File')
                            if start == False:
                                dict_df = {k: [] for k in properties}
                                start = True
                            for property in properties:
                                if property == 'File':
                                    dict_df[property].append(basename)
                                else:
                                    dict_df[property].append(final_dict[property])
                elif category == 'sp_ieea':
                     start = False
                     for file_name in dict.keys():
                            basename = self.file_base(file_name)
                            final_dict = dict[file_name]
                            ie = final_dict['ie']['E']
                            ea = final_dict['ea']['E']
                            if start == False:
                                dict_df = {k: [] for k in ['File', 'sp_ie','sp_ea']}
                                start=True
                            dict_df['File'].append(basename)
                            dict_df['sp_ie'].append(ie)
                            dict_df['sp_ea'].append(ea)
                elif category == 'ad_ieea':
                     start = False
                     for file_name in dict.keys():
                            basename = self.file_base(file_name)
                            final_dict = dict[file_name]
                            ie = final_dict['ie']['E']
                            ea = final_dict['ea']['E']
                            if start == False:
                                dict_df = {k: [] for k in ['File', 'ad_ie','ad_ea']}
                                start=True
                            dict_df['File'].append(basename)
                            dict_df['ad_ie'].append(ie)
                            dict_df['ad_ea'].append(ea)
                dict_df = pd.DataFrame(dict_df)
                if mol_df.empty:
                    mol_df = dict_df
                else:
                    mol_df = mol_df.merge(dict_df,how='left', on='File')
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