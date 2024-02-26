######################################################.
#        This file stores the IE & EA class           #
######################################################.


import sys, os
import time
import cclib as cc
from collections import defaultdict
from argument_parser import load_variables


class ie_ea:
    """
    Class containing all the functions for the  vertical ie and ea module related to Gaussian output files
    """

    def __init__(self, data, create_dat=True, **kwargs):

        start_time_overall = time.time()
        # load default and user-specified variables
        self.args = load_variables(kwargs, "IE_EA", create_dat=create_dat)
        self.data = data

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
                f"\nTime Collecting IE & EA data: {elapsed_time} seconds\n"
            )
            self.args.log.finalize()

    def get_data(self):
        mydict = lambda: defaultdict(mydict)
        file_data = mydict()

        for file_name in self.data.keys():
            ie_data, ea_data = None, None
            if "ie" in self.data[file_name].keys():
                ie_data = self.parse_cc_data(file_name, self.data[file_name]["ie"])
            if "ea" in self.data[file_name].keys():
                ea_data = self.parse_cc_data(file_name, self.data[file_name]["ea"])

            if ie_data != None and ea_data != None:
                self.args.log.write(
                    f"Reading information for IE & EA data from {file_name}\n"
                )
                file_data[file_name]["ie"]["E"] = ie_data.scfenergies[-1]
                file_data[file_name]["ea"]["E"] = ea_data.scfenergies[-1]
            else:
                self.args.log.write(
                    f"Skipping file {file_name} as one either IE or EA didnt exist\n"
                )

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
