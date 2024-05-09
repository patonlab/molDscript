######################################################.
#        This file stores the spc class               #
######################################################.


import sys, os
import time
import cclib as cc
from collections import defaultdict
from dftdescp.argument_parser import load_variables
import numpy as np
from openbabel import openbabel as ob
from rdkit import Chem

eV_to_hartree = 0.0367493


class spc:
    """
    Class containing all the functions for the opt module related to Gaussian output files
    """

    def __init__(self, data, create_dat=True, **kwargs):

        start_time_overall = time.time()
        # load default and user-specified variables
        self.args = load_variables(kwargs, "SPC", create_dat=create_dat)
        self.data = data

        if len(self.data.keys()) == 0:
            self.args.log.write(
                f"\nx  Could not find files to obtain information for single point correction"
            )
            self.args.log.finalize()
            sys.exit()
        else:
            self.file_data = self.get_data()

        if create_dat:
            elapsed_time = round(time.time() - start_time_overall, 2)
            self.args.log.write(f"   --- Single Point Energy Collection complete in {elapsed_time} seconds\n")
            self.args.log.finalize()

    def get_data(self):
        mydict = lambda: defaultdict(mydict)
        file_data = mydict()

        self.args.log.write(
                    f"   --- Single Point Energy Collection starting"
                )

        for file_name in self.data.keys():
            nickname = file_name

            spc_data = self.parse_cc_data(file_name, self.data[file_name])
            file_name = self.data[file_name]

            if self.args.program == 'gaussian':
                if list(self.data.keys()).index(nickname) == 0:
                    self.args.log.write(f"   Functional used: {spc_data.metadata['functional']}")
                    self.args.log.write(f"   Basis set used: {spc_data.metadata['basis_set']}")
            file_name = nickname
                
            
            self.args.log.write(f"o  Parsing SPC Energy Data from {os.path.basename(file_name)}")
            file_data[file_name]["spc_energy"] = (
                spc_data.scfenergies[-1] * eV_to_hartree)

        return file_data

    def parse_cc_data(self, file_name, file):

        ### parse data
        parser = cc.io.ccopen(file)
        try:
            cc_data = parser.parse()
        except:
            self.args.log.write(
                f"\nx  Could not parse {file_name} to obtain spc energy information"
            )
            cc_data = None

        return cc_data