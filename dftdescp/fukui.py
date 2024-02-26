######################################################.
#        This file stores the FUKUI class            #
######################################################.


import sys, os
import time
import cclib as cc
from collections import defaultdict
from argument_parser import load_variables


class fukui:
    """
    Class containing all the functions for the fukui module related to Gaussian output files
    """

    def __init__(self, data, create_dat=True, **kwargs):

        start_time_overall = time.time()
        # load default and user-specified variables
        self.args = load_variables(kwargs, "FUKUI", create_dat=create_dat)
        self.data = data

        if len(self.data.keys()) == 0:
            self.args.log.write(
                f"\nx  Could not find files to obtain information for calculating Fukui Coefficients"
            )
            self.args.log.finalize()
            sys.exit()
        else:
            self.file_data = self.get_data()

        if create_dat:
            elapsed_time = round(time.time() - start_time_overall, 2)
            self.args.log.write(
                f"\nTime Collecting FUKUI data: {elapsed_time} seconds\n"
            )
            self.args.log.finalize()

    def get_data(self):
        mydict = lambda: defaultdict(mydict)
        file_data = mydict()

        for file_name in self.data.keys():
            neutral_data, oxidized_data, reduced_data = None, None, None
            if "neutral" in self.data[file_name].keys():
                neutral_data = self.parse_cc_data(
                    file_name, self.data[file_name]["neutral"]
                )
            if "oxidized" in self.data[file_name].keys():
                oxidized_data = self.parse_cc_data(
                    file_name, self.data[file_name]["oxidized"]
                )
            if "reduced" in self.data[file_name].keys():
                reduced_data = self.parse_cc_data(
                    file_name, self.data[file_name]["reduced"]
                )

            if neutral_data != None and oxidized_data != None and reduced_data != None:
                self.args.log.write(
                    f"Reading information for FUKUI data from {file_name}\n"
                )
                file_data[file_name]["neutral"]["atomcharges"]["natural"] = (
                    neutral_data.atomcharges["natural"]
                )
                file_data[file_name]["neutral"]["atomcharges"]["cm5"] = (
                    neutral_data.atomcharges["cm5"]
                )
                file_data[file_name]["neutral"]["atomcharges"]["mulliken"] = (
                    neutral_data.atomcharges["mulliken"]
                )
                file_data[file_name]["neutral"]["atomcharges"]["hirsfeld"] = (
                    neutral_data.atomcharges["hirsfeld"]
                )

                file_data[file_name]["oxidized"]["atomcharges"]["natural"] = (
                    oxidized_data.atomcharges["natural"]
                )
                file_data[file_name]["oxidized"]["atomcharges"]["cm5"] = (
                    oxidized_data.atomcharges["cm5"]
                )
                file_data[file_name]["oxidized"]["atomcharges"]["mulliken"] = (
                    oxidized_data.atomcharges["mulliken"]
                )
                file_data[file_name]["oxidized"]["atomcharges"]["hirsfeld"] = (
                    oxidized_data.atomcharges["hirsfeld"]
                )

                file_data[file_name]["neutral"]["atomcharges"]["natural"] = (
                    reduced_data.atomcharges["natural"]
                )
                file_data[file_name]["neutral"]["atomcharges"]["cm5"] = (
                    reduced_data.atomcharges["cm5"]
                )
                file_data[file_name]["neutral"]["atomcharges"]["mulliken"] = (
                    reduced_data.atomcharges["mulliken"]
                )
                file_data[file_name]["neutral"]["atomcharges"]["hirsfeld"] = (
                    reduced_data.atomcharges["hirsfeld"]
                )
            else:
                self.args.log.write(
                    f"Skipping file {file_name} as one either neutral, oxidized or reduced didnt exist\n"
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

        try:  # npa
            start_npop = None
            lines = open(file, "r").readlines()
            for i, line in enumerate(lines):
                if line.find("-    Spin") > -1:
                    start_npop = i + 3
            if start_npop != None:
                nat_charges = []
                end_npop = start_npop + len(cc_data.atomnos)
                for i in range(start_npop, end_npop):
                    nat_charges.append(float(lines[i].split()[2]))
                cc_data.atomcharges["natural"] = nat_charges
        except:
            cc_data.atomcharges["natural"] = None

        try:  # hirsfeld & cm5
            start_npop = None
            lines = open(file, "r").readlines()
            for i, line in enumerate(lines):
                if (
                    line.find(
                        "Hirshfeld charges, spin densities, dipoles, and CM5 charges using IRadAn=      5:"
                    )
                    > -1
                ):
                    start_npop = i + 2
            if start_npop != None:
                hers_charges = []
                cm5_charges = []
                end_npop = start_npop + len(cc_data.atomnos)
                for i in range(start_npop, end_npop):
                    hers_charges.append(float(lines[i].split()[2]))
                    cm5_charges.append(float(lines[i].split()[-1]))
                cc_data.atomcharges["hirsfeld"] = hers_charges
                cc_data.atomcharges["cm5"] = cm5_charges
        except:
            cc_data.atomcharges["hirsfeld"] = None
            cc_data.atomcharges["cm5"] = None

        return cc_data
