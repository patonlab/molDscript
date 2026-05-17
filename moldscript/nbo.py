######################################################.
#        This file stores the NBO class               #
######################################################.


import sys, os
import time
import datetime
from moldscript.argument_parser import load_variables
from moldscript.utils import (
    format_timedelta,
    get_filename,
    initiate_data_dict,
    parse_npa_charges,
    progress_iter,
    record_cpu_time,
    safe_parse,
)
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
            self.data_dict = initiate_data_dict(self.data, logger=self.args.log)
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
        for i, file_name in progress_iter(self.data.keys(), self.args.log):
            nbo_data = self.parse_cc_data(file_name, self.data[file_name])

            if i == 1 and nbo_data is not None:
                self.args.log.write(f"   Package used: {nbo_data.metadata['package']} {nbo_data.metadata['package_version']}")
                try:
                    nbo_version = self.parse_nbo_version(self.data[file_name])
                    self.args.log.write(f"   NBO version used: {nbo_version}")
                    self.args.log.write(f"   Functional used: {nbo_data.metadata['functional']}")
                    self.args.log.write(f"   Basis set used: {nbo_data.metadata['basis_set']}\n")
                except (AttributeError, KeyError):
                    pass

            # When opt and nbo files have different suffix-strip levels (e.g.
            # ORCA QCALC + NBO), self.data and self.data_dict are keyed
            # differently — keep both keys around so neither lookup breaks.
            matched_key = get_filename(file_name, self.data_dict) if nbo_data is not None else file_name
            if nbo_data is not None:
                self.args.log.write_only(f"o  Parsing NBO data from {file_name}")
                self.data_dict[matched_key]['atom']["natural_charge"] = nbo_data.atomcharges["natural"]
                self.data_dict[matched_key]['atom']["bond_orders"] = nbo_data.bondorders
                if nbo_data.bondorders_matrix != []:
                    self.data_dict[matched_key]['bond']["bond_order_matrix"] = nbo_data.bondorders_matrix
            else:
                self.args.log.write(f"Skipping file {file_name} as NBO data didnt exist\n")

            cpu_times = nbo_data.metadata.get("cpu_time") if nbo_data and hasattr(nbo_data, "metadata") else None
            self.module_cpu_seconds += record_cpu_time(self.data_dict, matched_key, self.data[file_name], cpu_times)

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
    def parse_cc_data(self, file_name, file):
        cc_data = safe_parse(file, logger=self.args.log, context="NBO data")
        if cc_data is None:
            return None
        try:
            setattr(cc_data, "bondorders", self.bondorders(file, cc_data))
        except Exception:
            setattr(cc_data, "bondorders", None)
        setattr(cc_data, "bondorders_matrix", self.bondorders_matrix(file, cc_data))
        try:
            cc_data.atomcharges["natural"] = parse_npa_charges(file, len(cc_data.atomnos), which="last")
        except Exception:
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

