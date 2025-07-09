######################################################.
#          This file stores functions used           #
#                in multiple modules                 #
######################################################.

import os
import ast
from pathlib import Path
import glob
import datetime
import numpy as np
import cclib as cc
import moldscript.xyz2mol as xyz2mol
from rdkit import Chem

k_B_hartree = 3.1668114e-6  # hartree/K
J_TO_AU = 4.184 * 627.509541 * 1000.0  # UNIT CONVERSION
eV_to_hartree = 0.0367493

# class for logging
class Logger:
    """
    Class that wraps a file object to abstract the logging.
    """

    # Class Logger to write output to a file
    def __init__(self, filein, append, suffix="dat", verbose=True):
        if verbose:
            self.log = open(f"{filein}_{append}.{suffix}", "w")
        else:
            self.log = ''

    def write(self, message):
        """
        Appends a newline character to the message and writes it into the file.

        Parameters
        ----------
        message : str
           Text to be written in the log file.
        """
        try:
            self.log.write(f"{message}\n")
        except AttributeError:
            pass
        print(f"{message}")

    def write_only(self, message):
        """
        Appends a newline character to the message and writes it into the file, but does not print it.

        Parameters
        ----------
        message : str
           Text to be written in the log file.
        """
        try:
            self.log.write(f"{message}\n")
        except AttributeError:
            pass

    def finalize(self):
        """
        Closes the file
        """
        try:
            self.log.close()
        except AttributeError:
            pass

def add_cpu_times(file_data):
    ''' add cpu times for all files'''
    total_cpu = datetime.timedelta(0)
    for filename in file_data.keys():
        total_cpu += file_data[filename]['cpu_time']
    
    return total_cpu
def initiate_data_dict(data):
    """
    Initiates a data dictionary to store all the data from the files.
    """
    print(f"Initializing data parsing with SMILES and geometry data")
    data_dict = {}
    for i, file_name in enumerate(data.keys()):
        nickname = file_name
        data_dict[file_name] = dict()
        data_dict[file_name]["mol"] = dict()
        data_dict[file_name]["atom"] = dict()
        data_dict[file_name]["bond"] = dict()
        parsed_data = parse_cc_data(file_name, data[file_name])
        try:
            mol = xyz2mol.xyz2mol(parsed_data.atomnos.tolist(), parsed_data.atomcoords[-1].tolist(), charge=parsed_data.charge)[0]
            smi = Chem.MolToSmiles(mol)
        except:
            print("Encountered an issue with the mol embedding. Skipping smiles string.")
        data_dict[file_name]["mol"]["smiles"] = smi if 'smi' in locals() else ''
        data_dict[file_name]["atom"]["atomnos"] = parsed_data.atomnos
        data_dict[file_name]["bond"]["bond_length"] = parsed_data.bond_data_matrix
        data_dict[file_name]["mol"]["scfenergy"] = (parsed_data.scfenergies[-1] * eV_to_hartree)
        data_dict['CPU_time'] = []
    return data_dict

def format_lists(value):
    '''
    Transforms strings into a list
    '''

    if not isinstance(value, list):
        try:
            value = ast.literal_eval(value)
        except (SyntaxError, ValueError):
            # this line fixes issues when using "[X]" or ["X"] instead of "['X']" when using lists
            value = value.replace('[',']').replace(',',']').replace("'",']').split(']')
            while('' in value):
                value.remove('')
    return value
def bond_data_matrix(data):
        try:
            coords = data.atomcoords[-1]
        except:
            coords = data
        bond_data_matrix_list = []
        for atom1 in range(len(coords)):
            row = []
            for atom2 in range(len(coords)):
                p1 = np.array(coords[atom1])
                p2 = np.array(coords[atom2])
                squared_dist = np.sum((p1 - p2) ** 2, axis=0)
                dist = np.sqrt(squared_dist)
                row.append(dist)
            bond_data_matrix_list.append(row)
        return bond_data_matrix_list
def parse_cc_data(file_name, file):
        try:              
            parser = cc.io.ccopen(file)
            cc_data = parser.parse()
            setattr(cc_data, "bond_data_matrix", bond_data_matrix(cc_data))
        except:
            print(
                f"\nx  Could not parse {file_name}")
            raise SystemExit(f"Error parsing {file_name}. Ensure the file is a valid cclib file.")
        return cc_data
def get_files(value):

        if value[-1]=='/':
            value = value[:-1]
        if (Path(f"{value}").exists() and os.getcwd() not in f"{value}"):
            list_of_val_log = glob.glob(f"{os.getcwd()}/{value}/*.log")
            list_of_val_out = glob.glob(f"{os.getcwd()}/{value}/*.out")
        else:
            list_of_val_log = glob.glob(f"{value}/*.log")
            list_of_val_out = glob.glob(f"{value}/*.out")
        length_out = len(list_of_val_out)
        length_log = len(list_of_val_log)
        if length_log >= length_out:
            list_of_val = list_of_val_log
        else:
            list_of_val = list_of_val_out
        return list_of_val
 
def find_nth(haystack: str, needle: str, n: int) -> int:
    start = haystack.find(needle)
    while start >= 0 and n > 1:
        start = haystack.find(needle, start+len(needle))
        n -= 1
    return start
def get_filename(fullname, dd):
    flist = list(dd.keys())
    tempname = fullname
    for i in range(fullname.count("_")+1):
        try:
            findex = flist.index(tempname)
            keyname = flist[findex]
            return keyname
        except:
            tempname = tempname.rsplit("_", 1)[0]
            print(tempname)
    print(
        f"Error processing file {fullname}. Ensure consistent naming as described in the docs."
    )
    raise SystemExit


