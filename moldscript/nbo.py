######################################################.
#        This file stores the NBO class               #
######################################################.


import sys, os
import time
import datetime
import cclib as cc
from moldscript.argument_parser import load_variables
from moldscript.utils import (
    initiate_data_dict,
    record_cpu_time,
    format_timedelta,
    resolve_data_key,
    run_file_jobs,
    cpu_times_seconds,
)


def _read_output_lines(file):
    with open(file, "r") as outfile:
        return outfile.readlines()


def _parse_nbo_version_from_lines(lines):
    version = None
    for line in lines:
        if line.find("******* NBO") > -1:
            version = " ".join(line.split()[1:3])
    return version


def _nbo_bondorders_from_lines(lines):
    start_wiberg, end_wiberg = None, None
    for i, line in enumerate(lines):
        if line.find("Wiberg bond index, Totals by atom:") > -1:
            start_wiberg = i + 4
        if line.find("NBI: Natural Binding Index (NCU strength parameters)") > -1:
            end_wiberg = i - 2

    if start_wiberg is None or end_wiberg is None:
        return None
    wiberg_bos = []
    for i in range(start_wiberg, end_wiberg):
        wiberg_bos.append(float(lines[i].split()[2]))
    return wiberg_bos


def _nbo_bondorders_matrix_from_lines(lines, bondorders):
    if not bondorders:
        return []

    start_wiberg_ind, end_wiberg_ind = None, None
    for i, line in enumerate(lines):
        if line.find("Wiberg bond index matrix ") > -1:
            start_wiberg_ind = i + 2
        if line.find("Wiberg bond index") > -1:
            end_wiberg_ind = i - 1

    wiberg_bos_matrix = []
    if start_wiberg_ind is not None and end_wiberg_ind is not None:
        for i in range(start_wiberg_ind, end_wiberg_ind):
            if lines[i].find("Atom") > -1:
                for j, atom_idx in enumerate(lines[i].split()):
                    wbo_ind = []
                    if atom_idx != "Atom":
                        for k in range(i + 2, i + 2 + len(bondorders)):
                            wbo_ind.append(lines[k].split()[j + 1])
                        wiberg_bos_matrix.append(wbo_ind)
    return wiberg_bos_matrix


def _nbo_npa_from_lines(lines, atom_count):
    start_npop = None
    for i, line in enumerate(lines):
        if line.find(" Summary of Natural Population Analysis:") > -1:
            start_npop = i + 6

    if start_npop is None:
        return None
    nat_charges = []
    end_npop = start_npop + atom_count
    for i in range(start_npop, end_npop):
        nat_charges.append(float(lines[i].split()[2]))
    return nat_charges


def _parse_nbo_job(job):
    file_name, source_path, matched_name = job
    try:
        parser = cc.io.ccopen(source_path)
        nbo_data = parser.parse()
        lines = _read_output_lines(source_path)
        bondorders = _nbo_bondorders_from_lines(lines)
        bondorders_matrix = _nbo_bondorders_matrix_from_lines(lines, bondorders)
        natural_charge = _nbo_npa_from_lines(lines, len(nbo_data.atomnos))
        return {
            "file_name": file_name,
            "source_path": source_path,
            "matched_name": matched_name,
            "natural_charge": natural_charge,
            "bondorders": bondorders,
            "bondorders_matrix": bondorders_matrix,
            "nbo_version": _parse_nbo_version_from_lines(lines),
            "metadata": getattr(nbo_data, "metadata", {}),
            "cpu_times": nbo_data.metadata.get("cpu_time") if hasattr(nbo_data, "metadata") else None,
            "error": None,
        }
    except BaseException as exc:
        return {
            "file_name": file_name,
            "source_path": source_path,
            "matched_name": matched_name,
            "error": f"Could not parse {file_name} to obtain NBO information: {exc}",
        }


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
            self.data_dict = initiate_data_dict(
                self.data,
                logger=self.args.log,
                workers=self.args.workers,
            )
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
        jobs = []
        for file_name in self.data.keys():
            source_path = self.data[file_name]
            matched_name = resolve_data_key(file_name, self.data_dict, module_name="NBO", logger=self.args.log)
            jobs.append((file_name, source_path, matched_name))

        for idx, result in enumerate(
            run_file_jobs(jobs, _parse_nbo_job, workers=self.args.workers, logger=self.args.log)
        ):
            if result.get("error"):
                self.args.log.write(f"\nx  {result['error']}")
                raise SystemExit

            if idx == 0:
                metadata = result["metadata"]
                try:
                    self.args.log.write(f"   Package used: {metadata['package']} {metadata['package_version']}")
                    self.args.log.write(f"   NBO version used: {result['nbo_version']}")
                    self.args.log.write(f"   Functional used: {metadata['functional']}")
                    self.args.log.write(f"   Basis set used: {metadata['basis_set']}\n")
                except:
                    pass

            self.args.log.write_only(f"o  Parsing NBO data from {result['file_name']}")
            matched_name = result["matched_name"]
            self.data_dict[matched_name]['atom']["natural_charge"] = result["natural_charge"]
            self.data_dict[matched_name]['atom']["bond_orders"] = result["bondorders"]
            if result["bondorders_matrix"] != []:
                self.data_dict[matched_name]['bond']["bond_order_matrix"] = result["bondorders_matrix"]

            self.module_cpu_seconds += cpu_times_seconds(result["cpu_times"])
            record_cpu_time(self.data_dict, matched_name, result["source_path"], result["cpu_times"])

        return self.data_dict

    def parse_nbo_version(self, file):
        return _parse_nbo_version_from_lines(_read_output_lines(file))
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

        lines = _read_output_lines(file)
        setattr(cc_data, "bondorders", _nbo_bondorders_from_lines(lines))
        setattr(cc_data, "bondorders_matrix", _nbo_bondorders_matrix_from_lines(lines, cc_data.bondorders))
        cc_data.atomcharges["natural"] = _nbo_npa_from_lines(lines, len(cc_data.atomnos))
        return cc_data

    def bondorders(self, file, cc_data):
        return _nbo_bondorders_from_lines(_read_output_lines(file))

    def bondorders_matrix(self, file, cc_data):
        return _nbo_bondorders_matrix_from_lines(_read_output_lines(file), cc_data.bondorders)

    def npa_data(self, file, cc_data):
        return _nbo_npa_from_lines(_read_output_lines(file), len(cc_data.atomnos))


