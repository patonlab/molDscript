######################################################.
#        This file stores the min_max class            #
######################################################.

import pandas as pd
import numpy as np

class min_max:
    def __init__(self, cut=0.95, temp = 298.15, spc=False, energies =None, prefix=''):
        self.prefix = prefix
        self.threshold = 1 - cut
        self.temp = temp
        self.spc = spc
        self.weight_dict = {}
        self.energies = energies
        self.mol_min_max_range()
        self.atom_min_max_range()
        self.bond_min_max_range()


    def mol_min_max_range(self):
        ensemble_mol_csv = str(self.prefix) + 'min_max_range_molecule_level.csv'
        mol_df = pd.read_csv(str(self.prefix) + 'molecule_level.csv')
        mol_df = pd.merge(mol_df, self.energies, on='filename')
        
        print('\u25A1  CALCULATING MIN, MAX, AND RANGE FOR MOLECULE-LEVEL DESCRIPTORS INTO {}'.format(ensemble_mol_csv))
        basenames = mol_df['filename'].str.split('_conf').str[0].unique()
        result_df = pd.DataFrame()       
        for name in basenames:
            tempdf = mol_df[mol_df['filename'].str.contains(name)]
            boltzmann_constant_hartree = 3.1668114e-6  # Hartree per Kelvin
            tempdf['Exponent'] = -tempdf['scfenergy'] / (boltzmann_constant_hartree * self.temp)
            # Shift exponents for numerical stability
            tempdf['Exponent_shifted'] = tempdf['Exponent'] - tempdf['Exponent'].min()
            tempdf['Boltzmann_weight'] = np.exp(tempdf['Exponent_shifted'])
            # Normalize weights
            Z = tempdf['Boltzmann_weight'].sum()
            tempdf['Boltzmann_weight_normalized'] = tempdf['Boltzmann_weight'] / Z
            # Filter conformers based on threshold
            threshold = self.threshold  # Assume threshold is defined in your class
            tempdf = tempdf[tempdf['Boltzmann_weight_normalized'] >= threshold]
            if tempdf.empty:
                print(f"All conformers of {name} have Boltzmann weights below the threshold. Skipping.")
                continue  # Skip to next molecule
            # Update weight dictionary with conformers that passed the threshold
            filtered_weights = pd.Series(tempdf['Boltzmann_weight_normalized'].values, index=tempdf['filename']).to_dict()
            self.weight_dict.update(filtered_weights)
            # Get numerical and non-numerical columns
            numerical_cols = tempdf.select_dtypes(include=[np.number]).columns.tolist()
            non_numerical_cols = tempdf.select_dtypes(exclude=[np.number]).columns.tolist()
            # Exclude intermediate calculation columns
            intermediate_cols = ['Exponent', 'Exponent_shifted', 'Boltzmann_weight', 'Boltzmann_weight_normalized', 'scfenergy']
            numerical_cols = [col for col in numerical_cols if col not in intermediate_cols]
            # Prepare output_dict
            output_dict = {}
            if len(tempdf) == 1:
                # Only one conformer
                single_values = tempdf.iloc[0][numerical_cols]
                for col in numerical_cols:
                    output_dict[col + '_min'] = single_values[col]
                    output_dict[col + '_max'] = single_values[col]
                    output_dict[col + '_range'] = 0
            else:
                # Multiple conformers
                min_values = tempdf[numerical_cols].min()
                max_values = tempdf[numerical_cols].max()
                range_values = max_values - min_values
                for col in numerical_cols:
                    output_dict[col + '_min'] = min_values[col]
                    output_dict[col + '_max'] = max_values[col]
                    output_dict[col + '_range'] = range_values[col]
            # Keep the first value for non-numerical columns, excluding original numerical columns
            first_row = tempdf.iloc[0]
            for col in non_numerical_cols:
                if col in numerical_cols:
                    continue  # Exclude original numerical columns
                output_dict[col] = first_row[col]
            # Set filename to name
            output_dict['filename'] = name
            # Create output DataFrame
            output_df = pd.DataFrame([output_dict])
            result_df = pd.concat([result_df, output_df], ignore_index=True)
        columns_order = ['filename', 'smiles'] + [col for col in result_df.columns if col not in ['filename', 'smiles']]
        result_df = result_df[columns_order]
        result_df = result_df.round(4)
        result_df.to_csv(ensemble_mol_csv, index=False)

    def atom_min_max_range(self):
        atom_df = pd.read_csv(str(self.prefix) + 'atom_level.csv')
        ensemble_atom_csv = str(self.prefix) + 'min_max_range_atom_level.csv'
        print('\u25A1  CALCULATING MIN, MAX, AND RANGE FOR ATOM-LEVEL DESCRIPTORS INTO {}'.format(ensemble_atom_csv))
        
        # Map the weights to the atomic DataFrame based on 'filename'
        # Only include conformers that passed the threshold (weights exist in self.weight_dict)
        atom_df = atom_df[atom_df['filename'].isin(self.weight_dict.keys())].copy()
        atom_df['Weight'] = atom_df['filename'].map(self.weight_dict)
        # Extract the base molecule name
        atom_df['basename'] = atom_df['filename'].str.split('_conf').str[0]
        result_df = pd.DataFrame()
        for name in atom_df['basename'].unique():
            tempdf = atom_df[atom_df['basename'] == name]
            if tempdf.empty:
                print(f"No conformers of {name} passed the threshold. Skipping.")
                continue  # Skip to next molecule
            numerical_cols = tempdf.select_dtypes(include=[np.number]).columns.tolist()
            # Exclude columns that are not descriptors
            exclude_cols = ['Weight', 'atom_index']
            numerical_cols = [col for col in numerical_cols if col not in exclude_cols]
            grouped = tempdf.groupby('atom_index')
            atom_dicts = []  # Initialize a list to collect output_dicts
            for atom_idx, group in grouped:
                output_dict = {}
                if len(group) == 1:
                    # Only one conformer for this atom
                    single_values = group.iloc[0][numerical_cols]
                    for col in numerical_cols:
                        output_dict[col + '_min'] = single_values[col]
                        output_dict[col + '_max'] = single_values[col]
                        output_dict[col + '_range'] = 0
                else:
                    # Multiple conformers for this atom
                    min_values = group[numerical_cols].min()
                    max_values = group[numerical_cols].max()
                    range_values = max_values - min_values
                    for col in numerical_cols:
                        output_dict[col + '_min'] = min_values[col]
                        output_dict[col + '_max'] = max_values[col]
                        output_dict[col + '_range'] = range_values[col]
                # Keep non-numerical columns from the first occurrence
                non_numerical_cols = group.select_dtypes(exclude=[np.number]).columns.tolist()
                first_row = group.iloc[0]
                for col in non_numerical_cols:
                    if col in numerical_cols:
                        continue  # Exclude original numerical columns
                    output_dict[col] = first_row[col]
                # Add 'atom_index' back
                output_dict['atom_index'] = atom_idx
                # Append the output_dict to the list
                atom_dicts.append(output_dict)
            
            # Create a DataFrame from the list of dictionaries
            output_atomic_df = pd.DataFrame(atom_dicts)
            result_df = pd.concat([result_df, output_atomic_df], ignore_index=True)
        result_df['filename'] = result_df['basename']
        result_df = result_df.drop('basename', axis=1)
        columns_order = ['filename', 'atom_index', 'atom_type'] + [col for col in result_df.columns if col not in ['filename', 'atom_index', 'atom_type']]
        result_df = result_df[columns_order]
        result_df.to_csv(ensemble_atom_csv, index=False)

    def bond_min_max_range(self):
        bond_df = pd.read_csv(str(self.prefix) + 'bond_level.csv')
        ensemble_bond_csv =str(self.prefix) +  'min_max_range_bond_level.csv'
        print('\u25A1  CALCULATING MIN, MAX, AND RANGE FOR BOND-LEVEL DESCRIPTORS INTO {}'.format(ensemble_bond_csv))
        bond_df = bond_df[bond_df['filename'].isin(self.weight_dict.keys())].copy()
        bond_df['Weight'] = bond_df['filename'].map(self.weight_dict)
        bond_df['basename'] = bond_df['filename'].str.split('_conf').str[0]
        result_df = pd.DataFrame()
        for name in bond_df['basename'].unique():
            tempdf = bond_df[bond_df['basename'] == name]
            if tempdf.empty:
                print(f"No conformers of {name} passed the threshold. Skipping.")
                continue  # Skip to next molecule
            numerical_cols = tempdf.select_dtypes(include=[np.number]).columns.tolist()
            # Exclude columns that are not descriptors
            exclude_cols = ['Weight', 'atom1_idx', 'atom2_idx']
            numerical_cols = [col for col in numerical_cols if col not in exclude_cols]
            grouped = tempdf.groupby(['atom1_idx', 'atom2_idx'])
            bond_dicts = []
            for bond_indices, group in grouped:
                output_dict = {}
                if len(group) == 1:
                    # Only one conformer for this bond
                    single_values = group.iloc[0][numerical_cols]
                    for col in numerical_cols:
                        output_dict[col + '_min'] = single_values[col]
                        output_dict[col + '_max'] = single_values[col]
                        output_dict[col + '_range'] = 0
                else:
                    # Multiple conformers for this bond
                    min_values = group[numerical_cols].min()
                    max_values = group[numerical_cols].max()
                    range_values = max_values - min_values
                    for col in numerical_cols:
                        output_dict[col + '_min'] = min_values[col]
                        output_dict[col + '_max'] = max_values[col]
                        output_dict[col + '_range'] = range_values[col]
                # Keep non-numerical columns from the first occurrence
                non_numerical_cols = group.select_dtypes(exclude=[np.number]).columns.tolist()
                first_row = group.iloc[0]
                for col in non_numerical_cols:
                    if col in numerical_cols or col in ['Weight']:
                        continue  # Exclude original numerical columns and 'Weight'
                    output_dict[col] = first_row[col]
                # Add 'atom1_idx' and 'atom2_idx' back
                output_dict['atom1_idx'] = bond_indices[0]
                output_dict['atom2_idx'] = bond_indices[1]
                bond_dicts.append(output_dict)
            
            # Create a DataFrame from the list of dictionaries
            output_bond_df = pd.DataFrame(bond_dicts)
            # Add 'filename' as the molecule name
            output_bond_df['filename'] = name
            result_df = pd.concat([result_df, output_bond_df], ignore_index=True)
        
        result_df['filename'] = result_df['basename']
        result_df = result_df.drop('basename', axis=1)
        columns_order = ['filename', 'atom1_idx', 'atom1', 'atom2_idx', 'atom2'] + [col for col in result_df.columns if col not in ['filename', 'atom1_idx', 'atom1', 'atom2_idx', 'atom2']]
        result_df = result_df[columns_order]
        result_df.to_csv(ensemble_bond_csv, index=False)
