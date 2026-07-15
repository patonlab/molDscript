######################################################.
#        This file stores the get_df class            #
######################################################.
import pandas as pd
import numpy as np
pd.options.mode.chained_assignment = None
import math
import periodictable
import datetime
from moldscript.utils import append_run_log, emit


class get_df:
    """
    Class to create 3 .csv files containing paramaters
    """

    def __init__(self, data_dicts, substructure="", prefix="", bond_filter=False, no_mol=False, no_atom=False, no_bond=False, mol_vector=False):
        self.dd = data_dicts
        self.substructure = substructure
        self.prefix = prefix
        self.no_bond_filter = bond_filter
        if no_mol:
            emit('Skipping molecule-level descriptors', style="yellow")
        else:
            mol_df = self.get_mol_df()
        if no_bond:
            emit('Skipping bond-level descriptors', style="yellow")
        else:
            bond_df = self.get_bond_df()
        if no_atom:
            emit('Skipping atom-level descriptors', style="yellow")
        else:
            atom_df = self.get_atom_df()
        self.get_time()

    def get_mol_df(self):
        mol_csv = str(self.prefix) + "molecule_level.csv"
        emit("Writing molecule-level descriptor table to {}".format(mol_csv), style="cyan")
        filenames = list(self.dd.keys())
        filenames.remove("CPU_time")
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
            append_run_log(f"\t- {prop}")
        moldf.insert(0, "filename", col)
        moldf = moldf.round(4)
        self.energies = moldf[['filename', 'scfenergy']]
        moldf.drop('scfenergy', axis=1, inplace=True)
        moldf.to_csv(mol_csv, index=False)

    def get_bond_df(self):
        bond_csv = str(self.prefix) + "bond_level.csv"
        emit("Writing bond-level pair descriptor table to {}".format(bond_csv), style="cyan")
        filenames = list(self.dd.keys())
        filenames.remove("CPU_time")
        data = self.dd
        props = list(data[filenames[0]]["bond"].keys())
        for prop in props:
            append_run_log(f"\t- {prop}")
        bonddf = pd.DataFrame()
        for fname in filenames:
            atoms = np.array(data[fname]["atom"]["atomnos"])
            filedf = pd.DataFrame()
            for prop in props:
                matrix = np.array(data[fname]["bond"][prop])
                atom1_idx, atom2_idx = np.tril_indices_from(matrix, k=-1)
                values = matrix[atom1_idx, atom2_idx]
                fnames = np.array(fname for i in range(len(values)))
                num1 = atoms[atom1_idx.astype(int)]
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
                    filedf = pd.merge(filedf, tempdf, on=list(filedf.columns.intersection(tempdf.columns)))
                else:
                    filedf = tempdf
            if not bonddf.empty:
                bonddf = pd.concat([bonddf, filedf], axis=0)
            else:
                bonddf = filedf
        if self.substructure != '':
            emit(f'Filtering bond rows to user-defined substructure: {self.substructure}', style="cyan")
            final_df = pd.DataFrame()
            for filename in filenames:
                idx = self.dd[filename]['substructure']
                temp_df = bonddf.loc[bonddf['filename'] == filename]
                filtered_df = temp_df.loc[temp_df['atom1_idx'].isin(idx) | temp_df['atom2_idx'].isin(idx)]
                final_df = pd.concat([final_df, filtered_df])
            bonddf = final_df    
        bonddf = bonddf.round(4)
        if not self.no_bond_filter:
            if 'bond_order_matrix' in props:
                emit('Applying default bond filter: bond order >= 0.1', style="cyan")
                bonddf['bond_order_matrix'] = pd.to_numeric(bonddf['bond_order_matrix'], errors='coerce')
                bonddf = bonddf[bonddf['bond_order_matrix'] >= 0.1]
            else:
                emit('Applying default bond filter: bond length <= 3 angstroms', style="cyan")
                bonddf = bonddf[bonddf['bond_length'] <= 3]
        bonddf.to_csv(bond_csv, index=False)

    def get_atom_df(self):
        atom_csv = str(self.prefix) + "atom_level.csv"
        emit("Writing atom-level descriptor table to {}".format(atom_csv), style="cyan")
        filenames = list(self.dd.keys())
        filenames.remove("CPU_time")
        data = self.dd
        props = list(data[filenames[0]]["atom"].keys())
        calced_props = props.copy()
        calced_props.remove("atomnos")
        if calced_props != []:
            for prop in calced_props:
                append_run_log(f"\t- {prop}")
        atomdf = pd.DataFrame()
        for fname in filenames:
            atom_level_data = data[fname]["atom"]
            atoms = data[fname]["atom"]["atomnos"]
            tempdic = {}
            for prop in props:
                try:
                    values = atom_level_data[prop]
                    if values is None:
                        raise Exception
                except Exception:
                    values = [''] * len(atoms)
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
            emit(f'Filtering atom rows to user-defined substructure: {self.substructure}', style="cyan")
            final_df = pd.DataFrame()
            for filename in filenames:
                idx = self.dd[filename]['substructure']
                temp_df = atomdf.loc[atomdf['filename'] == filename]
                filtered_df = temp_df.loc[temp_df['atom_index'].isin(idx)]
                final_df = pd.concat([final_df, filtered_df])
            atomdf = final_df  
        try:
            filename = filenames[0]
            list(self.dd[filename]['sterics'].keys())
            emit(f'Merging steric buried-volume descriptors into {atom_csv}', style="cyan")
            steric_df = pd.DataFrame()
            for filename in filenames:
                props = list(self.dd[filename]['sterics'].keys())
                props.remove('atom_index')
                idxes = self.dd[filename]['sterics']['atom_index']
                fnames = [filename for i in idxes]
                tempdic = {}
                tempdic['filename'] = fnames
                tempdic['atom_index'] = idxes
                for prop in props:
                    values =  self.dd[filename]['sterics'][prop]
                    tempdic[str(prop)] = values
                temp_df = pd.DataFrame(tempdic)
                steric_df = pd.concat([steric_df, temp_df])
            for prop in props:
                    append_run_log(f"\t- {prop}")
            atomdf = pd.merge(atomdf, steric_df, on=['filename', 'atom_index'], how='outer')
        except:pass

        columns_order = ['filename', 'atom_index', 'atom_type'] + [col for col in atomdf.columns if col not in ['filename', 'atom_index', 'atom_type']]
        atomdf = atomdf[columns_order]
        atomdf = atomdf.map(lambda x: np.nan if isinstance(x, str) and x.strip() == '' else x)
        atomdf = atomdf.round(4)
        atomdf.to_csv(atom_csv, index=False)

    def get_time(self):
        cpu_time = datetime.timedelta(0)
        fnames = list(self.dd.keys())
        fnames.remove('CPU_time')
        for file in fnames:
            cpu_time += self.dd[file]["CPU_time"]
        total_seconds = cpu_time.total_seconds()
        hours = total_seconds / 3600
        emit(f'Total parsed-job CPU time: {hours:.0f} hours', style="bold green")

    def get_atom_lab(self, num):
        label = periodictable.elements[num]
        return label
