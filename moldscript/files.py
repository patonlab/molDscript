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

    def __init__(self, calc, path, data_dict, suffix, create_dat=False, **kwargs):

        # load default and user-specified variables
        self.args = load_variables(kwargs, "FILES", create_dat=create_dat)
        self.calc = calc
        self.path = path
        self.suffix = suffix
        self.warn_suffix = False
        if type(self.path) == str:
            self.files = get_files(self.path)
        elif type(self.path) == list:
            self.files = []
            for i in self.path:
                for j in get_files(i):
                    self.files.append(j)
        self.data_dict = data_dict

        # if self.calc == "link":
        #     self.file_data = self.get_link()
        if self.calc == "opt":
            self.file_data = self.get_opt_or_substructure()
        if self.calc == "spc":
            self.file_data = self.get_spc()
        if self.calc == "charges":
            self.file_data = self.get_chg()
        if self.calc == "fmo":
            self.file_data = self.get_fmo()
        if self.calc == "nmr":
            self.file_data = self.get_nmr()
        if self.calc == "nbo":
            self.file_data = self.get_nbo()
        if self.calc == "fukui":
            self.file_data = self.get_fukui()

        if self.calc == "ad_ie_ea":
            self.file_data = self.get_ad_ie_ea()
        if self.calc == "substructure":
            self.file_data = self.get_opt_or_substructure()

        if create_dat:
            self.args.log.finalize()

    def get_opt_or_substructure(self):
        file_data = defaultdict(dict)
        for file in self.files:
            key_name = self.get_filename(file, self.suffix)
            file_data[key_name] = file
        return file_data

    def get_spc(self):
        file_data = defaultdict(dict)
        for file in self.files:
            key_name = self.get_filename(file, self.suffix)
            file_data[key_name] = file
        return file_data
    def get_chg(self):
        file_data = defaultdict(dict)
        for file in self.files:
            key_name = self.get_filename(file, self.suffix)
            file_data[key_name] = file
        return file_data
    def get_fmo(self):
        file_data = defaultdict(dict)
        for file in self.files:
            key_name = self.get_filename(file, self.suffix)
            file_data[key_name] = file
        return file_data
    def get_nmr(self):
        file_data = defaultdict(dict)
        for file in self.files:
            key_name = self.get_filename(file, self.suffix)
            file_data[key_name] = file
        return file_data
    def get_nbo(self):
        file_data = defaultdict(dict)
        for file in self.files:
            key_name = self.get_filename(file, self.suffix)
            file_data[key_name] = file
        return file_data

    def get_fukui(self):
        file_data = defaultdict(dict)
        cwd = os.getcwd()
        neut = os.path.join(cwd, self.path[0])
        red = os.path.join(cwd, self.path[1])
        ox = os.path.join(cwd, self.path[2])
        for file in self.files:

            if neut in file:
                key_name = self.get_filename(file, self.suffix[0])
                file_data[key_name]["neutral"] = file
            if red in file:
                key_name = self.get_filename(file, self.suffix[1])
                file_data[key_name]["reduced"] = file
            if ox in file:
                key_name = self.get_filename(file, self.suffix[2])
                file_data[key_name]["oxidized"] = file
        return file_data

    # def get_link(self):
    #     file_data = defaultdict(dict)
    #     for file in self.files:
    #         file_data[file] = file
    #     return file_data

    def get_filename(self, fullname, suffix):

        try:
            fullname = fullname.split("/")[-1]
        except:
            pass
        try:
            fullname = fullname.split(".log")[0]
        except:
            pass
        try:
            fullname = fullname.split(".out")[0]
        except:
            pass
        if suffix != '':
            fullname = fullname.split("_" + suffix)[0]
        elif self.warn_suffix == False:
            print(f"Warning: no suffix provided for {self.calc}, using full filename")
            print("If this is not intentional, it will cause issues with matching filenames")
            self.warn_suffix = True
        return fullname