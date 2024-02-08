######################################################.
#        This file stores the FILES class            #
######################################################.


import sys, os
import time
from dftdescp.utils import (
    get_files,
)
from dftdescp.argument_parser import (load_variables)

class files:
    """
    Class containing all the functions from the files module related to Gaussian output files
    """

    def __init__(self, calc, path, create_dat=True, **kwargs):

        start_time_overall = time.time()
        # load default and user-specified variables
        self.args = load_variables(kwargs, "FILES", create_dat=create_dat)
        self.calc = calc
        self.path = path

        self.files = get_files(self.path)

        if len(self.files) == 0:
            self.args.log.write(f'\nx  No files were found! Make sure that 1) the PATH to the files is correct, 2) the PATH doesn\'t start with "/", and 3) you use quotation marks if you are using *')
            self.args.log.finalize()
            sys.exit()
        
        if self.calc == 'link':
            self.args.log.write('Reading all linked files')
            self.file_data = self.get_link_data()
        if self.calc == 'nmr':
            self.args.log.write('Reading NMR files')
            self.file_data = self.get_nmr()
        if self.calc == 'nbo':
            self.args.log.write('Reading NBO files')
            self.file_data = self.get_nbo()
        if self.calc == 'fukui':
            self.args.log.write('Reading FUKUI files')
            self.file_data = self.get_fukui()
        if self.calc == 'sp_ie_ea':
            self.args.log.write('Reading Vertical IE & EA files')
            self.file_data = self.get_sp_ie_ea()
        if self.calc == 'ad_ie_ea':
            self.args.log.write('Reading Adiabatic IE & EA files')
            self.file_data = self.get_ad_ie_ea()

    def get_nmr(self):
        file_data = {}
        for file in self.files:
            if self.args.suffix_nmr in os.path.basename(file):
                data = open(file).readlines()
                file_data[file] = data
        return file_data

    def get_nbo(self):
        file_data = {}
        for file in self.files:
            if self.args.suffix_nbo in os.path.basename(file):
                data = open(file).readlines()
                file_data[file] = data
        return file_data

    def get_fukui(self):
        file_data = {}
        for file in self.files:
            if self.args.suffix_fukui in os.path.basename(file):
                key_name = os.path.basename(file).split(f'_{self.args.suffix_fukui}')
                data = open(file).readlines()
                file_data[(key_name[0],file)] = data
            if self.args.suffix_fukui_n in os.path.basename(file):
                key_name = os.path.basename(file).split(f'_{self.args.suffix_fukui_n}')
                data = open(file).readlines()
                file_data[(key_name[0],file)] = data
            if self.args.suffix_fukui_p in os.path.basename(file):
                key_name = os.path.basename(file).split(f'_{self.args.suffix_fukui_p}')
                data = open(file).readlines()
                file_data[(key_name[0],file)] = data
        return file_data

    def get_sp_ie_ea(self):
        file_data = {}
        for file in self.files:
            if self.args.suffix_sp_ie in os.path.basename(file):
                key_name = os.path.basename(file).split(f'_{self.args.suffix_sp_ie}')
                data = open(file).readlines()
                file_data[(key_name[0],file)] = data
            if self.args.suffix_sp_ea in os.path.basename(file):
                key_name = os.path.basename(file).split(f'_{self.args.suffix_sp_ea}')
                data = open(file).readlines()
                file_data[(key_name[0],file)] = data
        return file_data

    def get_ad_ie_ea(self):
        file_data = {}
        for file in self.files:
            if self.args.suffix_ad_ie in os.path.basename(file):
                key_name = os.path.basename(file).split(f'_{self.args.suffix_ad_ie}')
                data = open(file).readlines()
                file_data[(key_name[0],file)] = data
            if self.args.suffix_ad_ea in os.path.basename(file):
                key_name = os.path.basename(file).split(f'_{self.args.suffix_ad_ea}')
                data = open(file).readlines()
                file_data[(key_name[0],file)] = data
        return file_data
    
    def get_link(self):
        file_data = {}
        for file in self.files:
                data = open(file).readlines()
                file_data[file] = data
        return file_data
        
