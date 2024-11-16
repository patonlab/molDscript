######################################################.
#        This file stores the spc class               #
######################################################.


import sys, os
import time
import cclib as cc
from collections import defaultdict
from moldscript.argument_parser import load_variables
from moldscript.utils import eV_to_hartree

class spc:
    """
    Class containing all the functions for the opt module related to Gaussian output files
    """

    def __init__(self, data, data_dict, create_dat=True,  **kwargs):
        
        start_time_overall = time.time()
        # load default and user-specified variables
        self.args = load_variables(kwargs, "SPC", create_dat=create_dat)
        self.data = data
        self.data_dict = data_dict
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

        self.args.log.write(
                    f"   --- Single Point Energy Collection starting"
                )

        for file_name in self.data.keys():
            if self.data[file_name].rsplit('.',1)[1] == 'log':
                self.spc_program = 'gaussian'
            elif self.data[file_name].rsplit('.', 1)[1] =='out':
                self.spc_program = 'orca'
            spc_data = self.parse_cc_data(file_name, self.data[file_name])
            
            full_filename = self.data[file_name]
            file_end = full_filename.split('/')[-1]
            filename = file_end.rsplit('_', 1)[0]

            if self.spc_program == 'gaussian':
                if list(self.data.keys()).index(file_name) == 0:
                    self.args.log.write(f"   Functional used: {spc_data.metadata['functional']}")
                    self.args.log.write(f"   Basis set used: {spc_data.metadata['basis_set']}")
            spc_data.scfenergies[-1] * eV_to_hartree
            self.args.log.write(f"o  Parsing SPC Energy Data from {os.path.basename(file_name)}")
            self.data_dict[filename]['mol']['spc_energy'] = (
                spc_data.scfenergies[-1] * eV_to_hartree)
        return self.data_dict

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