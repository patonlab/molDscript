######################################################.
#        This file stores the IE & EA class           #
######################################################.


import sys, os
import time
import cclib as cc
from collections import defaultdict
from moldscript.argument_parser import load_variables
from moldscript.utils import eV_to_hartree

class ie_ea:
    """
    Class containing all the functions for the  vertical ie and ea module related to Gaussian output files
    """

    def __init__(self, calc, data, data_dict, create_dat=True, **kwargs):

        start_time_overall = time.time()
        # load default and user-specified variables
        self.args = load_variables(kwargs, "IE_EA", create_dat=create_dat)
        self.data = data
        self.calc = calc 
        self.data_dict = data_dict

        if len(self.data.keys()) == 0:
            self.args.log.write(
                f"\nx  Could not find files to obtain information for calculating IE and EA"
            )
            self.args.log.finalize()
            sys.exit()
        else:
            self.file_data = self.get_data()

        if create_dat:
            elapsed_time = round(time.time() - start_time_overall, 2)
            self.args.log.write(
                f"   --- IE & EA Parameter Collection complete in {elapsed_time} seconds\n"
            )
            self.args.log.finalize()

    def get_data(self):
        mydict = lambda: defaultdict(mydict)
        file_data = mydict()
        first = False
        self.args.log.write(
                    f"   --- IE & EA Parameter Collection starting"
                )

        for file_name in self.data.keys():
            ie_data, ea_data = None, None

            if "oxidized" in self.data[file_name].keys():
                ie_data = self.parse_cc_data(file_name, self.data[file_name]["oxidized"])
                if first == False and self.args.program=='gaussian':
                    self.args.log.write(f"   Functional used: {ie_data.metadata['functional']}")
                    self.args.log.write(f"   Basis set used: {ie_data.metadata['basis_set']}")
                    first = True
            if "reduced" in self.data[file_name].keys():
                ea_data = self.parse_cc_data(file_name, self.data[file_name]["reduced"])
                if first == False and self.args.program=='gaussian':
                    self.args.log.write(f"   Functional used: {ea_data.metadata['functional']}")
                    self.args.log.write(f"   Basis set used: {ea_data.metadata['basis_set']}")
                    first = True
            if ie_data != None and ea_data != None:
                self.args.log.write(
                    f"o  Parsing IE & EA data from {file_name}"
                )
                self.data_dict[file_name]['mol']['oxidized_energy'] = ie_data.scfenergies[-1]*eV_to_hartree
                self.data_dict[file_name]['mol']['reduced_energy'] = ea_data.scfenergies[-1]*eV_to_hartree
                
            else:
                self.args.log.write(
                    f"x  Skipping file {file_name} as either IE or EA doest not exist!"
                )

        return self.data_dict

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
