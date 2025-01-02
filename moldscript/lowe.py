######################################################.
#  This file stores the lowest energy class    #
######################################################.

import pandas as pd
import numpy as np


class lowe:
    def __init__(self, temp = 298.15, spc=False, prefix='', energies=None ):

        self.prefix = prefix
        self.energies = energies
        self.low_confs = []
        self.mol_lowe()
        self.atom_lowe()
        self.bond_lowe()

    def mol_lowe(self):
        ensemble_mol_csv = str(self.prefix) + 'lowest_energy_molecule_level.csv'
        mol_df = pd.read_csv(str(self.prefix) + 'molecule_level.csv')
        mol_df = pd.merge(mol_df, self.energies, on='filename')
        print('\n')
        print('\u25A1  INCLUDING ONLY LOWEST ENERGY CONFORMER MOL DATA INTO {}'.format(ensemble_mol_csv))
        basenames = mol_df['filename'].str.split('_conf').str[0].unique()
        weighted_df = pd.DataFrame()

        for name in basenames:
            tempdf = mol_df[mol_df['filename'].str.contains(name)]
            mine = tempdf['scfenergy'].min()
            outdf = tempdf[tempdf['scfenergy'] == mine]

            self.low_confs.append(outdf['filename'].iloc[0])
            weighted_df = pd.concat([weighted_df, outdf])
 
        columns_order = ['filename', 'smiles'] + [col for col in weighted_df.columns if col not in ['filename', 'smiles']]
        weighted_df = weighted_df[columns_order]  
        weighted_df = weighted_df.drop('scfenergy', axis=1)  
        weighted_df = weighted_df.round(4)  
        weighted_df.to_csv(ensemble_mol_csv, index=False)

    def atom_lowe(self):
        atom_df = pd.read_csv(str(self.prefix) + 'atom_level.csv')  
        ensemble_atom_csv =str(self.prefix) +  'lowest_energy_atom_level.csv'
        print('\u25A1  INCLUDING ONLY LOWEST ENERGY CONFORMER ATOM DATA INTO {}'.format(ensemble_atom_csv))
        # Map the weights to the atomic DataFrame based on 'filename

        weighted_df = atom_df[atom_df['filename'].isin(self.low_confs)]

        weighted_df['filename'] = [k.rsplit('_conf',1)[0] for k in weighted_df['filename']]
        columns_order = ['filename', 'atom_index', 'atom_type'] + [col for col in weighted_df.columns if col not in ['filename', 'atom_index', 'atom_type']]
        weighted_df = weighted_df[columns_order]
        weighted_df = weighted_df.round(4)  
        weighted_df.to_csv(ensemble_atom_csv, index=False)


    def bond_lowe(self):
        bond_df = pd.read_csv(str(self.prefix) + 'bond_level.csv')
        ensemble_bond_csv =str(self.prefix) +  'lowest_energy_bond_level.csv'
        print('\u25A1  INCLUDING ONLY LOWEST ENERGY CONFORMER BOND DATA INTO {}\n'.format(ensemble_bond_csv))
        # Map the weights to the atomic DataFrame based on 'filename'
        weighted_df = bond_df[bond_df['filename'].isin(self.low_confs)]
        weighted_df['filename'] = [k.rsplit('_conf',1)[0] for k in weighted_df['filename']]
        columns_order = ['filename', 'atom1_idx', 'atom1', 'atom2_idx', 'atom2'] + [col for col in weighted_df.columns if col not in ['filename', 'atom1_idx', 'atom1', 'atom2_idx', 'atom2']]
        weighted_df = weighted_df[columns_order]
        weighted_df = weighted_df.round(4)  
        weighted_df.to_csv(ensemble_bond_csv, index=False)
       