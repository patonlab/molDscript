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
        except:
            self.rad = ast.literal_eval(self.rad)
        total = len(self.data)
        last_step = 0
        for i, fname in enumerate(self.data.keys(), start=1):
            file_name = self.data[fname]
            percent = int((i / total) * 100) if total else 100
            step = percent // 5
            if step > last_step:
                for s in range(last_step + 1, step + 1):
                    # progress printed to stdout and written to file
                    self.args.log.write(f"Progress: {s * 5}% ({i}/{total})")
                last_step = step
            # per-file messages go to the module .dat only
            self.args.log.write_only(f'o  Calculating steric parameters for {fname}')
            if vall:
                ccdata = cclib.io.ccread(file_name)
                atoms = range(1,len(ccdata.atomnos)+1)
            else:
                self.args.log.write_only(f"\t - Limiting volume calculation to the following indices: {self.dd[fname]['substructure']}")
                try:
                    atoms = self.dd[fname]['substructure'] 
                except:
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
        self.args.log.write_only(
            f"Error processing file {fullname}. Ensure consistent naming as described in the docs."
        )

            

