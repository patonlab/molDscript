######################################################.
#        This file stores the NBO class               #
######################################################.


import sys, os
import time
import datetime
import cclib as cc
from moldscript.argument_parser import load_variables
from moldscript.utils import get_filename, initiate_data_dict, record_cpu_time, format_timedelta
class nbo:
    """
    Class containing all the functions for the NBO module related to Gaussian output files
    """

    def __init__(self, data, data_dict: dict, create_dat=True, **kwargs) -> None:

        start_time_overall = time.time()
        # load default and user-specified variables
        self.args = load_variables(kwargs, "NBO", create_dat=create_dat)
        self.data = data
        self.data_dict = data_dict
        self.module_cpu_seconds = 0.0
        if self.data_dict == {}:
            self.data_dict = initiate_data_dict(self.data)
        self.fnames = self.data_dict.keys()

        if len(self.data.keys()) == 0:
            self.args.log.write(f"\nx  Could not find files to obtain information for calculating NBO")
            self.args.log.finalize()
            sys.exit()
        else:
            self.file_data = self.get_data()

        if create_dat:
            elapsed_time = round(time.time() - start_time_overall, 2)
            module_cpu_td = datetime.timedelta(seconds=self.module_cpu_seconds)
            if self.module_cpu_seconds:
                self.args.log.write(f"\n   NBO calculations CPU time: {format_timedelta(module_cpu_td)}")
            self.args.log.write(f"-- NBO Parameter Collection complete in {elapsed_time} seconds\n")
            self.args.log.finalize()

    def get_data(self):

        self.args.log.write(f"-- NBO Parameter Collection starting")
        self.module_cpu_seconds = 0.0
        for i, file_name in enumerate(self.data.keys()):
            nbo_data = self.parse_cc_data(file_name, self.data[file_name])

            if i == 0:
                self.args.log.write(f"   Package used: {nbo_data.metadata['package']} {nbo_data.metadata['package_version']}")
                try:
                    nbo_version = self.parse_nbo_version(self.data[file_name])
                    self.args.log.write(f"   NBO version used: {nbo_version}")
                    self.args.log.write(f"   Functional used: {nbo_data.metadata['functional']}")
                    self.args.log.write(f"   Basis set used: {nbo_data.metadata['basis_set']}\n")
                except:
                    pass

            if nbo_data != None:
                self.args.log.write(f"o  Parsing NBO data from {file_name}")
                file_name = get_filename(file_name, self.data_dict)
                self.data_dict[file_name]['atom']["natural_charge"] = nbo_data.atomcharges["natural"]
                self.data_dict[file_name]['atom']["bond_orders"] = nbo_data.bondorders
                if nbo_data.bondorders_matrix != []:
                    self.data_dict[file_name]['bond']["bond_order_matrix"] = nbo_data.bondorders_matrix


            else:
                self.args.log.write(f"Skipping file {file_name} as NBO data didnt exist\n")

            cpu_times = nbo_data.metadata.get("cpu_time") if nbo_data and hasattr(nbo_data, "metadata") else None
            self.module_cpu_seconds += record_cpu_time(self.data_dict, file_name, self.data[file_name], cpu_times)

        return self.data_dict

    def parse_nbo_version(self, file):
        start_version = None
        outfile = open(file, "r")
        lines = outfile.readlines()
        for i, line in enumerate(lines):
            if line.find("******* NBO") > -1:
                start_version = i

        if start_version != None:
            version = ' '.join(lines[start_version].split()[1:3])
        return version
    def get_filename(self):
        pass

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

        try:
            setattr(cc_data, "bondorders", self.bondorders(file, cc_data))
        except:
            setattr(cc_data, "bondorders", None)

        setattr(cc_data, "bondorders_matrix", self.bondorders_matrix(file, cc_data))

        try:
            cc_data.atomcharges["natural"] = self.npa_data(file, cc_data)
        except:
            cc_data.atomcharges["natural"] = None
        return cc_data

    def bondorders(self, file, cc_data):
        start_wiberg, end_wiberg = None, None
        outfile = open(file, "r")
        lines = outfile.readlines()
        for i, line in enumerate(lines):
            if line.find("Wiberg bond index, Totals by atom:") > -1:
                start_wiberg = i + 4
            if (
                line.find("NBI: Natural Binding Index (NCU strength parameters)")
                > -1
            ):
                end_wiberg = i - 2

        if start_wiberg != None and end_wiberg != None:
            wiberg_bos = []
            for i in range(start_wiberg, end_wiberg):
                wiberg_bos.append(float(lines[i].split()[2]))
        return wiberg_bos

    def bondorders_matrix(self, file, cc_data):

        start_wiberg_ind, end_wiberg_ind = None, None
        outfile = open(file, "r")
        lines = outfile.readlines()

        for i, line in enumerate(lines):
            if line.find("Wiberg bond index matrix ") > -1:
                start_wiberg_ind = i + 2
            if line.find("Wiberg bond index") > -1:
                end_wiberg_ind = i - 1
        wiberg_bos_matrix = []
        if start_wiberg_ind != None and end_wiberg_ind != None:

            for i in range(start_wiberg_ind, end_wiberg_ind):
                if lines[i].find("Atom") > -1:
                    for j, atom_idx in enumerate(lines[i].split()):
                        wbo_ind = []
                        if atom_idx != "Atom":
                            for k in range(i + 2, i + 2 + len(cc_data.bondorders)):
                                wbo_ind.append(lines[k].split()[j + 1])
                            wiberg_bos_matrix.append(wbo_ind)
        return wiberg_bos_matrix

    def npa_data(self, file, cc_data):
        start_npop = None
        outfile = open(file, "r")
        lines = outfile.readlines()
        for i, line in enumerate(lines):
            if line.find(" Summary of Natural Population Analysis:") > -1:
                start_npop = i + 6

        if start_npop != None:
            nat_charges = []
            end_npop = start_npop + len(cc_data.atomnos)
            for i in range(start_npop, end_npop):
                nat_charges.append(float(lines[i].split()[2]))
        return nat_charges


