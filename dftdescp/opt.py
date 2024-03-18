######################################################.
#        This file stores the OPT class               #
######################################################.


import sys, os
import time
import cclib as cc
from collections import defaultdict
from dftdescp.argument_parser import load_variables

eV_to_hartree = 0.0367493


class opt:
    """
    Class containing all the functions for the opt module related to Gaussian output files
    """

    def __init__(self, data, create_dat=True, **kwargs):

        start_time_overall = time.time()
        # load default and user-specified variables
        self.args = load_variables(kwargs, "OPT", create_dat=create_dat)
        self.data = data

        if len(self.data.keys()) == 0:
            self.args.log.write(
                f"\nx  Could not find files to obtain information optimization"
            )
            self.args.log.finalize()
            sys.exit()
        else:
            self.file_data = self.get_data()

        if create_dat:
            elapsed_time = round(time.time() - start_time_overall, 2)
            self.args.log.write(f"   --- Optimization Parameter Collection complete in {elapsed_time} seconds\n")
            self.args.log.finalize()

    def get_data(self):
        mydict = lambda: defaultdict(mydict)
        file_data = mydict()

        self.args.log.write(
                    f"   --- Optimization Parameter Collection starting"
                )

        for file_name in self.data.keys():
            opt_data = self.parse_cc_data(file_name, self.data[file_name])
            if list(self.data.keys()).index(file_name) == 0:
                self.args.log.write(f'Functional used: {opt_data.metadata['functional']}')
                self.args.log.write(f'Basis set used: {opt_data.metadata['basis_set']}')
            
            self.args.log.write(f"o  Parsing Energy & Thermochemistry Data from {os.path.basename(file_name)}")
            file_data[file_name]["opt"]["scfenergy"] = (
                opt_data.scfenergies[-1] * eV_to_hartree
            )
            file_data[file_name]["opt"]["enthalpy"] = opt_data.enthalpy
            file_data[file_name]["opt"]["freeenergy"] = opt_data.freeenergy
        return file_data

    def parse_cc_data(self, file_name, file):

        ### parse data
        parser = cc.io.ccopen(file)
        try:
            cc_data = parser.parse()
        except:
            self.args.log.write(
                f"\nx  Could not parse {file_name} to obtain information for calculating Fukui Coefficients"
            )
            cc_data = None

        return cc_data
