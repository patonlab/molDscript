######################################################.
#        This file stores the NMR class               #
######################################################.


import sys, os
import time
import cclib as cc
from collections import defaultdict
from moldscript.argument_parser import load_variables

class nmr:
    """
    Class containing all the functions for the NMR module related to Gaussian output files
    """

    def __init__(self, data, create_dat=True, **kwargs):

        start_time_overall = time.time()
        # load default and user-specified variables
        self.args = load_variables(kwargs, "NMR", create_dat=create_dat)
        self.data = data

        if len(self.data.keys()) == 0:
            self.args.log.write(
                f"\nx  Could not find files to obtain information for calculating NMR"
            )
            self.args.log.finalize()
            sys.exit()
        else:
            self.file_data = self.get_data()

        if create_dat:
            elapsed_time = round(time.time() - start_time_overall, 2)
            self.args.log.write(f"   --- NMR Parameter Collection complete in {elapsed_time} seconds\n")
            self.args.log.finalize()

    def get_data(self):
        mydict = lambda: defaultdict(mydict)
        file_data = mydict()

        self.args.log.write(
                    f"   --- NMR Parameter Collection starting"
                )
        for file_name in self.data.keys():
            try:
                nmr_data = self.parse_cc_data(file_name, self.data[file_name])
            except:
                nmr_data = None

            if list(self.data.keys()).index(file_name) == 0 and self.args.program=='gaussian':
                self.args.log.write(f"   Functional used: {nmr_data.metadata['functional']}")
                self.args.log.write(f"   Basis set used: {nmr_data.metadata['basis_set']}")
            if nmr_data != None:
                self.args.log.write(
                    f"o  Parsing NMR Shielding Tensors from {file_name}"
                )
                file_data[file_name]["nmr_shielding"] = nmr_data.nmr_shielding
                file_data[file_name]['atomnos'] = nmr_data.atomnos
            else:
                self.args.log.write(
                    f"!  Skipping {file_name} as NMR data not found"
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

        if self.args.program=='gaussian':
            try: setattr(cc_data, "nmr_shielding", self.gaussian_nmr_shielding(file))
            except: setattr(cc_data, "nmr_shielding", None)
        
        if self.args.program=='orca':
            try: setattr(cc_data, "nmr_shielding", self.orca_nmr_shielding(file, cc_data.natom))
            except: setattr(cc_data, "nmr_shielding", None)

        return cc_data

    def gaussian_nmr_shielding(self, file):
        outfile = open(file, "r")
        lines = outfile.readlines()
        for i in range(0, len(lines)):
            if lines[i].find("shielding tensors") > -1:
                start = i + 2
            if lines[i].find("End of Minotr F.D.") > -1:
                end = i - 1
        nmr_shielding = []
        for j in range(start, end - 1, 5):
            nmr = lines[j].split()[4]
            nmr_shielding.append(nmr)
        return nmr_shielding

    def orca_nmr_shielding(self, file, natoms):
        outfile = open(file, "r")
        lines = outfile.readlines()
        for i in range(0, len(lines)):
            if lines[i].find("CHEMICAL SHIELDING SUMMARY (ppm)") > -1:
                start = i + 6
        
        end = start + natoms
        nmr_shielding = []
        for j in range(start, end):
            nmr = lines[j].split()[2]
            nmr_shielding.append(nmr)
        return nmr_shielding
        