######################################################.
#        This file stores the IE & EA class           #
######################################################.


import sys, os
import time
import datetime
import cclib as cc
from collections import defaultdict
from moldscript.argument_parser import load_variables
from moldscript.utils import eV_to_hartree, add_cpu_times

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
            try:
                total_cpu = add_cpu_times(self.file_data)
                self.args.log.write(f"\n   IE & EA calculations complete in {total_cpu} seconds")
            except: pass
  
            self.args.log.write(
                f"-- IE & EA Parameter Collection complete in {elapsed_time} seconds\n"
            )
            self.args.log.finalize()

    def get_data(self):
        mydict = lambda: defaultdict(mydict)
        file_data = mydict()

        for i, file_name in enumerate(self.data.keys()):
            ie_data, ea_data = None, None

            if "ie" in self.data[file_name].keys():
                ie_data = self.parse_cc_data(file_name, self.data[file_name]["ie"])

            if "ea" in self.data[file_name].keys():
                ea_data = self.parse_cc_data(file_name, self.data[file_name]["ea"])
            
            if i == 0:
                rel_dir = self.data[file_name]["ie"].split(os.getcwd()+'/')[1].split(file_name)[0]    
                self.args.log.write(
                    f"-- IE & EA Parameter Collection from {rel_dir}"
                )
                self.args.log.write(f"   Package used: {ea_data.metadata['package']} {ea_data.metadata['package_version']}")
                self.args.log.write(f"   Functional used: {ea_data.metadata['functional']}")
                self.args.log.write(f"   Basis set used: {ea_data.metadata['basis_set']}\n")
            
            if ie_data != None and ea_data != None:
                self.args.log.write(
                    f"o  Parsing IE & EA data from {file_name}"
                )
                file_data[file_name]["ox"]["E"] = ie_data.scfenergies[-1] *eV_to_hartree
                file_data[file_name]["red"]["E"] = ea_data.scfenergies[-1] *eV_to_hartree
                
            else:
                self.args.log.write(
                    f"x  Skipping file {file_name} as either IE or EA doest not exist!"
                )

            file_data[file_name]['cpu_time'] = datetime.timedelta(0) # initialize cpu time
            for time in ie_data.metadata['cpu_time']:
                file_data[file_name]['cpu_time'] += time # add cpu time from IE
            for time in ea_data.metadata['cpu_time']:
                file_data[file_name]['cpu_time'] += time # add cpu time from EA

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
