######################################################.
#        This file stores the NMR class               #
######################################################.


import sys, os
import time
import datetime
from collections import defaultdict
from moldscript.argument_parser import load_variables
from moldscript.utils import (
    format_timedelta,
    get_filename,
    initiate_data_dict,
    progress_iter,
    record_cpu_time,
    safe_parse,
)


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
            self.data_dict = initiate_data_dict(self.data, logger=self.args.log)
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

        self.args.log.write(f"-- NMR Parameter Collection starting")
        self.module_cpu_seconds = 0.0
        for i, file_name in progress_iter(self.data.keys(), self.args.log):
            # Resolve the matched data_dict key separately from the file-side
            # key — when opt and nmr files have different suffix-strip levels
            # (e.g. ORCA QCALC + nmr) the two diverge, and self.data is keyed
            # by the file-side basename, not the matched dict key.
            matched_key = get_filename(file_name, self.data_dict, logger=self.args.log)
            nmr_data = self.parse_cc_data(file_name, self.data[file_name])
            if i == 1 and nmr_data is not None:
                try:
                    self.args.log.write(f"   Package used: {nmr_data.metadata['package']} {nmr_data.metadata['package_version']}")
                    self.args.log.write(f"   Functional used: {nmr_data.metadata['functional']}")
                    self.args.log.write(f"   Basis set used: {nmr_data.metadata['basis_set']}\n")
                except (AttributeError, KeyError):
                    pass
            if nmr_data is not None:
                self.args.log.write_only(f"o  Parsing NMR Shielding Tensors from {matched_key}")
                self.data_dict[matched_key]["atom"]["nmr_shielding"] = nmr_data.nmr_shielding
            else:
                self.args.log.write(f"!  Skipping {matched_key} as NMR data not found")

            cpu_times = nmr_data.metadata.get("cpu_time") if nmr_data and hasattr(nmr_data, "metadata") else None
            self.module_cpu_seconds += record_cpu_time(self.data_dict, matched_key, self.data[file_name], cpu_times)
        module_cpu_td = datetime.timedelta(seconds=self.module_cpu_seconds)
        if self.module_cpu_seconds:
            self.args.log.write(f"-- NMR CPU time: {format_timedelta(module_cpu_td)}")
        return self.data_dict

    def parse_cc_data(self, file_name, file):
        # Thin wrapper: shared cclib parse + per-package shielding-tensor scrape.
        cc_data = safe_parse(file, logger=self.args.log, context="NMR shielding tensors")
        if cc_data is None:
            return None
        package = cc_data.metadata['package'].lower()
        if package == "gaussian":
            cc_data.nmr_shielding = self.gaussian_nmr_shielding(file)
        elif package == "orca":
            cc_data.nmr_shielding = self.orca_nmr_shielding(file)
        else:
            raise ValueError(f"Unsupported NMR tensor program: {cc_data.metadata['package']}")
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


