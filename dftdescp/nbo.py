######################################################.
#        This file stores the NBO class               #
######################################################.


import sys, os
import time
import cclib as cc
from collections import defaultdict
from dftdescp.argument_parser import load_variables


class nbo:
    """
    Class containing all the functions for the NBO module related to Gaussian output files
    """

    def __init__(self, data, create_dat=True, **kwargs):

        start_time_overall = time.time()
        # load default and user-specified variables
        self.args = load_variables(kwargs, "NBO", create_dat=create_dat)
        self.data = data

        if len(self.data.keys()) == 0:
            self.args.log.write(
                f"\nx  Could not find files to obtain information for calculating NBO"
            )
            self.args.log.finalize()
            sys.exit()
        else:
            self.file_data = self.get_data()

        if create_dat:
            elapsed_time = round(time.time() - start_time_overall, 2)
            self.args.log.write(f"   --- NBO Parameter Collection complete in {elapsed_time} seconds\n")
            self.args.log.finalize()

    def get_data(self):
        mydict = lambda: defaultdict(mydict)
        file_data = mydict()

        self.args.log.write(
                    f"   --- NBO Parameter Collection starting"
                )

        for file_name in self.data.keys():
            try:
                nbo_data = self.parse_cc_data(file_name, self.data[file_name])
            except:
                nbo_data = None
            if list(self.data.keys()).index(file_name) == 0:
                self.args.log.write(f"   Functional used: {nbo_data.metadata['functional']}")
                self.args.log.write(f"   Basis set used: {nbo_data.metadata['basis_set']}")
            if nbo_data != None:
                self.args.log.write(
                    f"o  Parsing NBO & NPA data from {file_name}"
                )
                file_data[file_name]["charges"]["npa"] = nbo_data.atomcharges["natural"]
                file_data[file_name]["bond_orders"] = nbo_data.bondorders
                file_data[file_name]["bond_order_matrix"] = nbo_data.bondorders_matrix
            else:
                self.args.log.write(
                    f"Skipping file {file_name} as NBO data didnt exist\n"
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

        try:  # grab Wiberg Bond Orders summed for each atom
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
                setattr(cc_data, "bondorders", wiberg_bos)
            else:
                setattr(cc_data, "bondorders", None)
        except:
            setattr(cc_data, "bondorders", None)

        try:  # grab individual Wiberg Bond Orders for each atom
            start_wiberg_ind, end_wiberg_ind = None, None
            outfile = open(file, "r")
            lines = outfile.readlines()
            for i, line in enumerate(lines):
                if line.find("Wiberg bond index matrix ") > -1:
                    start_wiberg_ind = i + 2
                if line.find("Wiberg bond index") > -1:
                    end_wiberg_ind = i - 1
            if start_wiberg_ind != None and end_wiberg_ind != None:
                wiberg_bos_matrix = []
                for i in range(start_wiberg_ind, end_wiberg_ind):
                    if lines[i].find("Atom") > -1:
                        for j, atom_idx in enumerate(lines[i].split()):
                            wbo_ind = []
                            if atom_idx != "Atom":
                                for k in range(i + 2, i + 2 + len(cc_data.bondorders)):
                                    wbo_ind.append(lines[k].split()[j + 1])
                                wiberg_bos_matrix.append(wbo_ind)

                setattr(cc_data, "bondorders_matrix", wiberg_bos_matrix)
            else:
                setattr(cc_data, "bondorders_matrix", None)
        except:
            setattr(cc_data, "bondorders_matrix", None)

        try:
            start_npop = None
            outfile = open(file, "r")
            lines = outfile.readlines()
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

        return cc_data
