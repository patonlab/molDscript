######################################################.
#        This file stores the sterimol descriptor class      #
######################################################.


import dbstep.Dbstep as db
from rdkit import Chem
import os
import cclib
import pandas as pd
from openbabel import openbabel as ob
import time
import ast

class sterics:
    def __init__(self, opt_data, data_dict, volume, vall, radii=3):
        t1 = time.time()
        self.data = opt_data    
        self.dd = data_dict
        self.rad = radii
        if vall != False:
            self.get_params(vall=vall)
        else:
            self.get_params(volume=volume)
        elapsed_time = round(time.time() - t1, 2)
        print(f"-- Steric Parameter Collection complete in {elapsed_time} seconds\n")
    def get_params(self, vall=False, volume=False):
        try:
            self.rad = [float(self.rad)]
        except:
            self.rad = ast.literal_eval(self.rad)

        for i, fname in enumerate(self.data.keys()):
            file_name = self.data[fname]
            print(f'o  Calculating steric parameters for {fname}')
            if vall:
                ccdata = cclib.io.ccread(file_name)
                atoms = range(1,len(ccdata.atomnos)+1)
            else:
                print(f"\t - Limiting volume calculation to the following indices: {self.dd[fname]['substructure']}")
                try:
                    atoms = self.dd[fname]['substructure'] 
                except:
                    print('\tPlease include a --substructure or specify --vall to get all buried volumes.')
                    break           
            self.dd[fname]['sterics'] = {}              
            for radius in self.rad:
                radius = float(radius)
                self.get_vbur(radius, atoms, fname, file_name)
            


    def get_vbur(self, radius, indexes, fname, file_name):
        
        bv_list = []
        idx_list = []
        for i in indexes:
            mol = db.dbstep(file_name, atom1=i, b=True, quiet=True, r=radius)
            bv = mol.bur_vol
            bv_list.append(bv)
            idx_list.append(i)
        self.dd[fname]['sterics'][f'buried_volume_r_{radius}'] = bv_list   
        self.dd[fname]['sterics'][f'atom_index'] = idx_list
        #remove the temp files created by dbstep
        f_woe = file_name.rsplit('.', 1)[0] 
        tmp_py = f_woe+ '_steric.py'
        tmp_xyz = f_woe + '_transform.xyz'
        os.remove(tmp_py)
        os.remove(tmp_xyz)

    def get_filename(self, fullname):
        flist = list(self.dd.keys())
        try:
            fullname = fullname.split("/")[-1]
        except:
            pass
        tempname = fullname
        for i in range(fullname.count("_")):
            try:
                findex = flist.index(tempname)
                keyname = flist[findex]
                return keyname
            except:
                tempname = tempname.rsplit("_", 1)[0]
        print(
            f"Error processing file {fullname}. Ensure consistent naming as described in the docs."
        )

            

