######################################################.
#        This file stores the sterimol descriptor class      #
######################################################.


import dbstep.Dbstep as db
from rdkit import Chem
import os
import cclib
import pandas as pd
import time
import ast
from moldscript.argument_parser import load_variables
from moldscript.utils import progress_iter

class sterics:
    def __init__(self, opt_data, data_dict, volume, vall, radii=3):
        t1 = time.time()
        # create a module logger so messages go to MOLDSCRIPT_STERICS.dat
        self.args = load_variables({}, "STERICS", create_dat=True)
        self.data = opt_data
        self.dd = data_dict
        self.rad = radii
        if vall != False:
            self.get_params(vall=vall)
        else:
            self.get_params(volume=volume)
        elapsed_time = round(time.time() - t1, 2)
        self.args.log.write_only(f"-- Steric Parameter Collection complete in {elapsed_time} seconds\n")
    def get_params(self, vall=False, volume=False):
        try:
            self.rad = [float(self.rad)]
        except (TypeError, ValueError):
            self.rad = ast.literal_eval(self.rad)
        for i, fname in progress_iter(self.data.keys(), self.args.log):
            file_name = self.data[fname]
            self.args.log.write_only(f'o  Calculating steric parameters for {fname}')
            if vall:
                ccdata = cclib.io.ccread(file_name)
                atoms = range(1,len(ccdata.atomnos)+1)
            else:
                self.args.log.write_only(f"\t - Limiting volume calculation to the following indices: {self.dd[fname]['substructure']}")
                try:
                    atoms = self.dd[fname]['substructure']
                except KeyError:
                    # Note: the log.write_only call above already accesses
                    # ['substructure'] unconditionally and will KeyError first
                    # if the key is missing, so this branch is effectively dead.
                    self.args.log.write_only('\tPlease include a --substructure or specify --vall to get all buried volumes.')
                    break
            self.dd[fname]['sterics'] = {}              
            for radius in self.rad:
                radius = float(radius)
                self.args.log.write(f'\t - Calculating buried volume for {fname} at radius {radius}.\n\tThis may take a moment...')
                self.get_vbur(radius, atoms, fname, file_name)
            


    def get_vbur(self, radius, indexes, fname, file_name):
        
        bv_list = []
        idx_list = []
        for i in indexes:
            mol = db.dbstep(file_name, atom1=i, b=True, quiet=True, r=radius, grid=0.1, noH=True, scalevdw=1.17, exclude=str(i))
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

