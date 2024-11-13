######################################################.
#        This file stores the sterimol descriptor class      #
######################################################.


import dbstep.Dbstep as db
from rdkit import Chem
import os
import cclib
import pandas as pd

class sterics:
    def __init__(self, file='../tests/gaussian/QCALC/success/Ac2_rdkit_conf_1.log', vbur=True, smol=False):
        self.file = file
        self.vbur = vbur
        self.smol = smol
        if self.vbur:
            self.volumes, self.numbers= self.get_vbur()

        
    def get_vbur(self):
        ccdata = cclib.io.ccread(self.file)
        indexes = ccdata.atomnos
        id_list = []
        bv_list = []
        type_list = []
        for i in range(len(indexes)):
            index = i + 1
            mol = db.dbstep(self.file, atom1=index, b=True, quiet=True)
            bv = mol.bur_vol
            type_list.append(ccdata.atomnos[i])
            bv_list.append(bv)
        fname = os.path.splitext(self.file)[0]
        tmp_py = fname + '_steric.py'
        tmp_xyz = fname + '_transform.xyz'
        os.remove(tmp_py)
        os.remove(tmp_xyz)
        return (bv_list, type_list) 
        
            

