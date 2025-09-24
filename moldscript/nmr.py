######################################################.
#        This file stores the NMR class               #
######################################################.


import sys, os
import time
import datetime
import cclib as cc
from collections import defaultdict
from moldscript.argument_parser import load_variables
from moldscript.utils import initiate_data_dict, record_cpu_time, format_timedelta


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
        self.module_cpu_seconds = 0.0
        if self.data_dict == {}:
            self.data_dict = initiate_data_dict(self.data)
        self.flist = list(data_dicts.keys())

        if len(self.data.keys()) == 0:
            self.args.log.write(f"\nx  Could not find files to obtain information for calculating NMR")
            self.args.log.finalize()
            sys.exit()
        else:
            self.file_data = self.get_data()

        if create_dat:
            elapsed_time = round(time.time() - start_time_overall, 2)
            self.args.log.write(f"-- NMR Parameter Collection complete in {elapsed_time} seconds\n")
            self.args.log.finalize()

    def get_data(self):
        mydict = lambda: defaultdict(mydict)
        file_data = mydict()

        print(f"-- NMR Parameter Collection starting")
        self.module_cpu_seconds = 0.0
        for i, file_name in enumerate(self.data.keys()):
            # try:
            file_name = self.get_filename(file_name)
            nmr_data = self.parse_cc_data(file_name, self.data[file_name])
            # except:
            #     nmr_data = None
            try:
                if i == 0:
                    print(f"   Package used: {nmr_data.metadata['package']} {nmr_data.metadata['package_version']}")
                    print(f"   Functional used: {nmr_data.metadata['functional']}")
                    print(f"   Basis set used: {nmr_data.metadata['basis_set']}\n")
            except:
                pass
            if nmr_data != None:
                print(f"o  Parsing NMR Shielding Tensors from {file_name}")
                self.data_dict[file_name]["atom"]["nmr_shielding"] = nmr_data.nmr_shielding
            else:
                print(f"!  Skipping {file_name} as NMR data not found")

            cpu_times = nmr_data.metadata.get("cpu_time") if nmr_data and hasattr(nmr_data, "metadata") else None
            self.module_cpu_seconds += record_cpu_time(self.data_dict, file_name, self.data[file_name], cpu_times)
        module_cpu_td = datetime.timedelta(seconds=self.module_cpu_seconds)
        if self.module_cpu_seconds:
            print(f"-- NMR CPU time: {format_timedelta(module_cpu_td)}")
        return self.data_dict

    def parse_cc_data(self, file_name, file):

        parser = cc.io.ccopen(file)
        cc_data = parser.parse()
        if cc_data.metadata['package'].lower() == "gaussian":
            cc_data.nmr_shielding = self.gaussian_nmr_shielding(file)
        elif cc_data.metadata['package'].lower() == "orca":
            cc_data.nmr_shielding = self.orca_nmr_shielding(file)
        else:
            raise ValueError(f"Unsupported nmr tensor program: {cc_data.metadata['package']}")


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
            if nmr == "Anisotropy": #Bad formatting in Gaussian
                nmr_line = lines[j].split()[3]
                nmr = nmr_line.strip("=")
            nmr = float(nmr)
            nmr_shielding.append(nmr)
        return nmr_shielding

    def orca_nmr_shielding(self, file):
        outfile = open(file, "r")
        lines = outfile.readlines()
        for i in range(0, len(lines)):
            if lines[i].find("CHEMICAL SHIELDING SUMMARY (ppm)") > -1:
                idx = i + 6
        nmr_shielding = []
        line = True
        for i in range(idx, len(lines)):
            line = lines[idx]
            if line.split() == []:
                break
            nmr = float(line.split()[2])
            nmr_shielding.append(nmr)
            idx += 1
        return nmr_shielding

    def get_filename(self, fullname):
        flist = list(self.data_dict.keys())
        tempname = fullname
        for i in range(fullname.count("_")+1):
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

