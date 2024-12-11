######################################################.
#        This file stores the NMR class               #
######################################################.


import sys, os
import time
import datetime
import cclib as cc
from collections import defaultdict
from moldscript.argument_parser import load_variables
from moldscript.utils import add_cpu_times


class nmr:
    """
    Class containing all the functions for the NMR module related to Gaussian output files
    """

    def __init__(self, data, data_dicts, create_dat=True, **kwargs):

        start_time_overall = time.time()
        # load default and user-specified variables
        self.args = load_variables(kwargs, "NMR", create_dat=create_dat)
        self.data = data
        self.data_dict = data_dicts
        self.flist = list(data_dicts.keys())

        if len(self.data.keys()) == 0:
            self.args.log.write(f"\nx  Could not find files to obtain information for calculating NMR")
            self.args.log.finalize()
            sys.exit()
        else:
            self.file_data = self.get_data()

        if create_dat:
            elapsed_time = round(time.time() - start_time_overall, 2)
            try:
                total_cpu = add_cpu_times(self.file_data)
                self.args.log.write(f"\n   NMR calculations complete in {total_cpu} seconds")
            except:
                pass
            self.args.log.write(f"-- NMR Parameter Collection complete in {elapsed_time} seconds\n")
            self.args.log.finalize()

    def get_data(self):
        mydict = lambda: defaultdict(mydict)
        file_data = mydict()

        self.args.log.write(f"-- NMR Parameter Collection starting")
        for i, file_name in enumerate(self.data.keys()):
            try:
                file_name = self.get_filename(file_name)
                nmr_data = self.parse_cc_data(file_name, self.data[file_name])
            except:
                nmr_data = None
            try:
                if i == 0:
                    self.args.log.write(f"   Package used: {nmr_data.metadata['package']} {nmr_data.metadata['package_version']}")
                    self.args.log.write(f"   Functional used: {nmr_data.metadata['functional']}")
                    self.args.log.write(f"   Basis set used: {nmr_data.metadata['basis_set']}\n")
            except:
                pass
            if nmr_data != None:
                self.args.log.write(f"o  Parsing NMR Shielding Tensors from {file_name}")
                self.data_dict[file_name]["atom"]["nmr_shielding"] = nmr_data.nmr_shielding
            else:
                self.args.log.write(f"!  Skipping {file_name} as NMR data not found")

            if self.data[file_name] in self.data_dict['CPU_time']:
                pass
            else:
                try: 
                    for time in nmr_data.metadata["cpu_time"]:
                        self.data_dict[file_name]["CPU_time"] += time  # add cpu time
                    self.data_dict["CPU_time"].append(self.data[file_name])
                except:
                    self.data_dict[file_name]["CPU_time"] = datetime.timedelta(0)  # initialize cpu time
                    for time in nmr_data.metadata["cpu_time"]:
                        self.data_dict[file_name]["CPU_time"] += time  # add cpu time
                    self.data_dict['CPU_time'].append(self.data[file_name])
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
        if file.rsplit(".", 1)[1] == "log":
            try:
                setattr(cc_data, "nmr_shielding", self.gaussian_nmr_shielding(file))
            except:
                setattr(cc_data, "nmr_shielding", None)

        elif file.rsplit(".", 1)[1] == "out":
            try:
                setattr(
                    cc_data,
                    "nmr_shielding",
                    self.orca_nmr_shielding(file, cc_data.natom),
                )
            except:
                setattr(cc_data, "nmr_shielding", None)

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
            nmr = float(lines[j].split()[4])
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
            nmr = float(lines[j].split()[2])
            nmr_shielding.append(nmr)
        return nmr_shielding

    def get_filename(self, fullname):
        flist = list(self.data_dict.keys())
        tempname = fullname
        for i in range(fullname.count("_")):
            try:
                findex = flist.index(tempname)
                keyname = flist[findex]
                return keyname
            except:
                tempname = tempname.rsplit("_", 1)[0]
                print(tempname)
        print(
            f"Error processing file {fullname}. Ensure consistent naming as described in the docs."
        )
        raise SystemExit
