######################################################.
#        This file stores the FILES class            #
######################################################.


import sys, os
import time
from moldscript.utils import (
    get_files,
)
from collections import defaultdict
from moldscript.argument_parser import load_variables


class files:
    """
    Class containing all the functions from the files module related to output files
    """

    def __init__(self, calc, path, create_dat=False, **kwargs):

        start_time_overall = time.time()
        # load default and user-specified variables
        self.args = load_variables(kwargs, "FILES", create_dat=create_dat)
        self.calc = calc
        self.path = path
        self.program = self.args.program
        self.files = get_files(self.path, self.program)

        if self.calc == "link":
            self.file_data = self.get_link()
        if self.calc == "opt":
            self.file_data = self.get_opt_or_substurcture()
        if self.calc == "spc":
            self.file_data = self.get_spc()
        if self.calc == "nmr":
            self.file_data = self.get_nmr()
        if self.calc == "nbo":
            self.file_data = self.get_nbo()
        if self.calc == "fukui":
            self.file_data = self.get_fukui()
        if self.calc == "sp_ie_ea":
            self.file_data = self.get_sp_ie_ea()
        if self.calc == "ad_ie_ea":
            self.file_data = self.get_ad_ie_ea()
        if self.calc == "substructure":
            self.file_data = self.get_opt_or_substurcture()

        if create_dat:
            self.args.log.finalize()

    def get_opt_or_substurcture(self):
        file_data = defaultdict(dict)
        for file in self.files:
            
            if self.args.suffix_opt == None:
                ftype = '.log'
                if self.args.program == 'orca':
                    ftype = '.out'
                key_name = os.path.basename(file).split(ftype)

                file_data[key_name[0]] = file
            elif self.args.suffix_opt in os.path.basename(file):
                key_name = os.path.basename(file).split(f"_{self.args.suffix_opt}")
                file_data[key_name[0]] = file
        return file_data
    
    def get_spc(self):
        file_data = defaultdict(dict)
        for file in self.files:            
            if self.args.suffix_spc in os.path.basename(file):
                key_name = os.path.basename(file).split(f"_{self.args.suffix_spc}")
                file_data[key_name[0]] = file
        return file_data

    def get_nmr(self):
        file_data = defaultdict(dict)
        for file in self.files:         
            if self.args.suffix_nmr in os.path.basename(file):
                key_name = os.path.basename(file).split(f"_{self.args.suffix_nmr}")
                file_data[key_name[0]] = file
        return file_data

    def get_nbo(self):
        file_data = defaultdict(dict)
        for file in self.files:
            if self.args.suffix_nbo in os.path.basename(file):
                key_name = os.path.basename(file).split(f"_{self.args.suffix_nbo}")
                file_data[key_name[0]] = file
        return file_data

    def get_fukui(self):
        file_data = defaultdict(dict)
        for file in self.files:
            if self.args.suffix_fukui in os.path.basename(file):
                key_name = os.path.basename(file).split(f"_{self.args.suffix_fukui}")
                file_data[key_name[0]]["neutral"] = file
            if self.args.suffix_fred in os.path.basename(file):
                key_name = os.path.basename(file).split(
                    f"_{self.args.suffix_fred}"
                )
                file_data[key_name[0]]["reduced"] = file
            if self.args.suffix_fox in os.path.basename(file):
                key_name = os.path.basename(file).split(f"_{self.args.suffix_fox}")
                file_data[key_name[0]]["oxidized"] = file
        return file_data

    def get_sp_ie_ea(self):
        file_data = defaultdict(dict)
        for file in self.files:
            if self.args.suffix_sp_ie in os.path.basename(file):
                key_name = os.path.basename(file).split(f"_{self.args.suffix_sp_ie}")
                file_data[key_name[0]]["ie"] = file
            elif self.args.suffix_sp_ea in os.path.basename(file):
                key_name = os.path.basename(file).split(f"_{self.args.suffix_sp_ea}")
                file_data[key_name[0]]["ea"] = file
        return file_data

    def get_ad_ie_ea(self):
        file_data = defaultdict(dict)
        for file in self.files:
            if self.args.suffix_ad_ie in os.path.basename(file):
                key_name = os.path.basename(file).split(f"_{self.args.suffix_ad_ie}")
                file_data[key_name[0]]["ie"] = file
            if self.args.suffix_ad_ea in os.path.basename(file):
                key_name = os.path.basename(file).split(f"_{self.args.suffix_ad_ea}")
                file_data[key_name[0]]["ea"] = file
        return file_data

    def get_link(self):
        file_data = defaultdict(dict)
        for file in self.files:
            file_data[file] = file
        return file_data
