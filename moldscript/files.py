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

    def __init__(self, calc, path, data_dict, create_dat=False, **kwargs):

        start_time_overall = time.time()
        # load default and user-specified variables
        self.args = load_variables(kwargs, "FILES", create_dat=create_dat)
        self.calc = calc
        self.path = path
        self.program = self.args.program
        if type(self.path) == str:
            self.files = get_files(self.path, self.program)
        elif type(self.path) == list:
            self.files = []
            for i in self.path:
                for j in get_files(i, self.program):
                    self.files.append(j)
        self.data_dict = data_dict

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
            ftype = ".log"
            if self.args.program == "orca":
                ftype = ".out"
            key_name = os.path.basename(file).split(ftype)
            file_data[key_name[0]] = file
        return file_data

    def get_spc(self):
        file_data = defaultdict(dict)
        for file in self.files:
            key_name = self.get_filename(file)
            file_data[key_name] = file
        return file_data

    def get_nmr(self):
        file_data = defaultdict(dict)
        for file in self.files:
            key_name = self.get_filename(file)
            file_data[key_name] = file
        return file_data

    def get_nbo(self):
        file_data = defaultdict(dict)
        for file in self.files:
            key_name = self.get_filename(file)
            file_data[key_name] = file
        return file_data

    def get_fukui(self):
        file_data = defaultdict(dict)
        neut = self.path[0]
        red = self.path[1]
        ox = self.path[2]
        for file in self.files:
            key_name = self.get_filename(file)

            if neut in file:
                file_data[key_name]["neutral"] = file
            if red in file:
                file_data[key_name]["reduced"] = file
            if ox in file:
                file_data[key_name]["oxidized"] = file
        return file_data

    def get_ad_ie_ea(self):
        file_data = defaultdict(dict)
        red = self.path[0]
        ox = self.path[1]
        for file in self.files:
            key_name = self.get_filename(file)
            if red in file:
                file_data[key_name]["reduced"] = file
            if ox in file:
                file_data[key_name]["oxidized"] = file
        return file_data

    def get_link(self):
        file_data = defaultdict(dict)
        for file in self.files:
            file_data[file] = file
        return file_data

    def get_filename(self, fullname):
        flist = list(self.data_dict.keys())
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
        print(
            f"Error processing file {fullname}. Ensure consistent naming as described in the docs."
        )
