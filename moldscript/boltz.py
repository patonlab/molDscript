######################################################.
#  This file stores the boltzmann weighting class    #
######################################################.

import pandas as pd
import numpy as np
boltzmann_constant = 3.1668114e-6

class boltz:
    def __init__(self, temp = 298.15,  prefix='', energies=None ):

        self.temp = temp
        self.prefix = prefix
        self.energies = energies
        self.weight_dict = {}
        self.mol_boltz()
        self.atom_boltz()
        self.bond_boltz()

    def mol_boltz(self):
        ensemble_mol_csv = str(self.prefix) + 'ensemble_molecule_level.csv'
        mol_df = pd.read_csv(str(self.prefix) + 'molecule_level.csv')
        mol_df = pd.merge(mol_df, self.energies, on='filename')

        print('\u25A1  AVERAGING MOLECULE-LEVEL DESCRIPTORS OVER CONFORMERS INTO {}'.format(ensemble_mol_csv))
        basenames = mol_df['filename'].str.split('_conf').str[0].unique()
        weighted_df = pd.DataFrame()
        for name in basenames:
            tempdf = mol_df[mol_df['filename'].str.split('_conf').str[0] == name]
            fullname = list(tempdf['filename'])[0]
            if len(tempdf) == 1:
                tempdf['filename'] = [name]
                output_df = tempdf
                self.weight_dict[fullname] = 1
            else:
                # Calculate Boltzmann weights
                # Make scfenergy values relative to the lowest one
                # work on a copy to avoid SettingWithCopy issues
                tempdf = tempdf.copy()
                tempdf['scfenergy'] = tempdf['scfenergy'] - tempdf['scfenergy'].min()
                tempdf['Exponent'] = -tempdf['scfenergy'] / (boltzmann_constant * self.temp)

                tempdf['Boltzmann_weight'] = np.exp(tempdf['Exponent'])
                # Normalize Boltzmann weights so they sum to 1
                tempdf['Boltzmann_weight'] = tempdf['Boltzmann_weight'] / tempdf['Boltzmann_weight'].sum()
                # print(name)
                # print(tempdf[['filename', 'scfenergy', 'Exponent', 'Boltzmann_weight']])
                # Normalize weights
                Z = tempdf['Boltzmann_weight'].sum()
                tempdf['Boltzmann_weight_normalized'] = tempdf['Boltzmann_weight'] / Z
                # Create a dictionary of weights with 'filename' as keys
                new_weights = pd.Series(tempdf['Boltzmann_weight_normalized'].values, index=tempdf['filename']).to_dict()
                self.weight_dict.update(new_weights)
                # Prepare the output row
                numerical_cols = tempdf.select_dtypes(include=[np.number]).columns.tolist()
                non_numerical_cols = tempdf.select_dtypes(exclude=[np.number]).columns.tolist()
                # Exclude intermediate calculation columns from numerical columns
                intermediate_cols = ['Exponent', 'Exponent_shifted', 'Boltzmann_weight', 'Boltzmann_weight_normalized']
                numerical_cols = [col for col in numerical_cols if col not in intermediate_cols]
                # Weight numerical columns
                weighted_values = {}
                for col in numerical_cols:
                    weighted_values[col] = (tempdf[col] * tempdf['Boltzmann_weight_normalized']).sum()
                # Keep the first value for non-numerical columns
                first_row = tempdf.iloc[0]
                for col in non_numerical_cols:
                    weighted_values[col] = first_row[col]
                # Create the output row as a DataFrame
                output_df = pd.DataFrame([weighted_values])
                output_df['filename'] = [name]
            weighted_df = pd.concat([weighted_df, output_df])
        columns_order = ['filename', 'smiles'] + [col for col in weighted_df.columns if col not in ['filename', 'smiles']]
        # Save self.weight_dict as a CSV file
        weights_df = pd.DataFrame(list(self.weight_dict.items()), columns=['filename', 'boltzmann_weight'])
        weights_df.to_csv(str(self.prefix) + 'boltzmann_weights.csv', index=False)
        weighted_df = weighted_df[columns_order]  
        weighted_df = weighted_df.drop('scfenergy', axis=1)  
        weighted_df = weighted_df.round(4)  
        weighted_df.to_csv(ensemble_mol_csv, index=False)

    def atom_boltz(self):
        try:
            atom_df = pd.read_csv(str(self.prefix) + 'atom_level.csv')
        except FileNotFoundError:
            print(f"atom_level.csv not found at {str(self.prefix) + 'atom_level.csv'}, skipping atom_boltz.")
            return
        
        ensemble_atom_csv =str(self.prefix) +  'ensemble_atom_level.csv'
        print('\u25A1  AVERAGING ATOM-LEVEL DESCRIPTORS OVER CONFORMERS INTO {}'.format(ensemble_atom_csv))
        # Map the weights to the atomic DataFrame based on 'filename'
        atom_df['Weight'] = atom_df['filename'].map(self.weight_dict)
        atom_df = atom_df.dropna(subset=['Weight'])
        atom_df['basename'] = atom_df['filename'].str.split('_conf').str[0]
        weighted_df = pd.DataFrame()
        for name in atom_df['basename'].unique():
            # select rows that exactly match this basename
            tempdf = atom_df[atom_df['basename'] == name]
            numerical_cols = tempdf.select_dtypes(include=[np.number]).columns.tolist()
            # remove columns if present (safe removal)
            for col_to_remove in ['Weight', 'atom_index']:
                if col_to_remove in numerical_cols:
                    numerical_cols.remove(col_to_remove)
            # Group by atom index to compute weighted averages
            grouped = tempdf.groupby('atom_index')
            weighted_avgs = grouped.apply(
    lambda x: pd.Series({col: (x[col] * x['Weight']).sum() / x['Weight'].sum() for col in numerical_cols}),
    include_groups=False)
            weighted_avgs.reset_index(inplace=True, drop=False)
            # For non-numerical columns, retain the first value in each group
            non_numerical_cols = tempdf.select_dtypes(exclude=[np.number]).columns.tolist()
            # Get the first occurrence of non-numerical columns per group
            first_values = grouped[non_numerical_cols].first().reset_index()

            # Merge weighted averages with first values
            output_atomic_df = pd.merge(weighted_avgs, first_values, on='atom_index')
            weighted_df = pd.concat([weighted_df, output_atomic_df])
        weighted_df['filename'] = weighted_df['basename']
        weighted_df = weighted_df.drop('basename', axis=1)
        columns_order = ['filename', 'atom_index', 'atom_type'] + [col for col in weighted_df.columns if col not in ['filename', 'atom_index', 'atom_type']]
        weighted_df = weighted_df[columns_order]
        weighted_df = weighted_df.round(4)  
        weighted_df.to_csv(ensemble_atom_csv, index=False)


    def bond_boltz(self):
        bond_df = pd.read_csv(str(self.prefix) + 'bond_level.csv')
        try:
            bond_df = pd.read_csv(str(self.prefix) + 'bond_level.csv')
        except FileNotFoundError:
            print(f"bond_level.csv not found at {str(self.prefix) + 'bond_level.csv'}, skipping bond_boltz.")
            return
        ensemble_bond_csv =str(self.prefix) +  'ensemble_bond_level.csv'
        print('\u25A1  AVERAGING BOND-LEVEL DESCRIPTORS OVER CONFORMERS INTO {}\n'.format(ensemble_bond_csv))
        # Map the weights to the atomic DataFrame based on 'filename'
        bond_df['Weight'] = bond_df['filename'].map(self.weight_dict)
        bond_df = bond_df.dropna(subset=['Weight'])
        bond_df['basename'] = bond_df['filename'].str.split('_conf').str[0]
        weighted_df = pd.DataFrame()
        for name in bond_df['basename'].unique():
            # select rows that exactly match this basename
            tempdf = bond_df[bond_df['basename'] == name]
            numerical_cols = tempdf.select_dtypes(include=[np.number]).columns.tolist()
            # safe removal of index/weight columns
            for col_to_remove in ['Weight', 'atom1_idx', 'atom2_idx']:
                if col_to_remove in numerical_cols:
                    numerical_cols.remove(col_to_remove)
            grouped = tempdf.groupby(['atom1_idx', 'atom2_idx'])
            weighted_avgs = grouped.apply(
    lambda x: pd.Series({col: (x[col] * x['Weight']).sum() / x['Weight'].sum() for col in numerical_cols}),
    include_groups=False)
            weighted_avgs.reset_index(inplace=True, drop=False)
            # For non-numerical columns, retain the first value in each group
            non_numerical_cols = tempdf.select_dtypes(exclude=[np.number]).columns.tolist()
            # Get the first occurrence of non-numerical columns per group
            first_values = grouped[non_numerical_cols].first().reset_index()
            # Merge weighted averages with first values
            output_bond_df = pd.merge(weighted_avgs, first_values, on=['atom1_idx', 'atom2_idx'])
            weighted_df = pd.concat([weighted_df, output_bond_df])
        weighted_df['filename'] = weighted_df['basename']
        weighted_df = weighted_df.drop('basename', axis=1)
        columns_order = ['filename', 'atom1_idx', 'atom1', 'atom2_idx', 'atom2'] + [col for col in weighted_df.columns if col not in ['filename', 'atom1_idx', 'atom1', 'atom2_idx', 'atom2']]
        weighted_df = weighted_df[columns_order]
        weighted_df = weighted_df.round(4)  
        weighted_df.to_csv(ensemble_bond_csv, index=False)
       