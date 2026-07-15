######################################################.
#        This file stores the spc class               #
######################################################.


import sys, os
import time
import datetime
import numpy as np
import cclib as cc
from moldscript.argument_parser import load_variables
from moldscript.utils import initiate_data_dict, record_cpu_time, format_timedelta, resolve_data_key

class charges:
    """
    Class containing all the functions for the charges module related to Gaussian output files
    """

    def __init__(self, data, data_dict, create_dat=True,  **kwargs):

        start_time_overall = time.time()
        # load default and user-specified variables
        self.args = load_variables(kwargs, "CHARGES", create_dat=create_dat)
        self.data = data
        self.data_dict = data_dict
        self.module_cpu_seconds = 0.0
        if self.data_dict == {}:
            self.data_dict = initiate_data_dict(self.data, logger=self.args.log)
        if len(self.data.keys()) == 0:
            self.args.log.write(f"\nx  Could not find files to obtain information for charge data")
            sys.exit()
        else:
            self.file_data = self.get_data()

        if create_dat:
            elapsed_time = round(time.time() - start_time_overall, 2)
            self.args.log.write(f"-- Charges Collection complete in {elapsed_time} seconds")
            self.args.log.finalize()

    def get_data(self):

        self.args.log.write(f"-- Charges Collection starting")
        self.module_cpu_seconds = 0.0
        total = len(self.data)
        last_step = 0
        for idx, file_name in enumerate(self.data.keys(), start=1):
            source_path = self.data[file_name]
            percent = int((idx / total) * 100) if total else 100
            step = percent // 5
            if step > last_step:
                for s in range(last_step + 1, step + 1):
                    self.args.log.write(f"Progress: {s * 5}% ({idx}/{total})")
                last_step = step
            chg_data = self.parse_cc_data(file_name, source_path)
            filename = self.get_filename(file_name)

            try:
                if list(self.data.keys()).index(file_name) == 0:
                    self.args.log.write(f"   Functional used: {chg_data.metadata['functional']}")
                    self.args.log.write(f"   Basis set used: {chg_data.metadata['basis_set']}")
            except:
                pass
            self.args.log.write_only(f"o  Parsing Charge Data from {os.path.basename(file_name)}")
            if len(chg_data.atomcharges.keys()) == 1 and 'mulliken' in chg_data.atomcharges:
                self.data_dict[filename]['atom']['mulliken_charge'] = chg_data.atomcharges['mulliken']
            else:
                for i in chg_data.atomcharges.keys():
                    if 'mulliken' not in i and 'sum' not in i:
                        self.data_dict[filename]['atom'][str(i)+'_charge'] = chg_data.atomcharges[i]

            for spin_type, spins in self.get_atom_spins(chg_data, source_path).items():
                self.data_dict[filename]['atom'][str(spin_type)+'_spin'] = spins

            cpu_times = chg_data.metadata.get("cpu_time") if chg_data and hasattr(chg_data, "metadata") else None
            self.module_cpu_seconds += record_cpu_time(self.data_dict, filename, source_path, cpu_times)
        module_cpu_td = datetime.timedelta(seconds=self.module_cpu_seconds)
        if self.module_cpu_seconds:
            self.args.log.write(f"-- Charges CPU time: {format_timedelta(module_cpu_td)}")
        return self.data_dict

    def parse_cc_data(self, file_name, file):

        parser = cc.io.ccopen(file)

        try:
            cc_data = parser.parse()
        except:
            self.args.log.write_only(
                f"\nx  Could not parse {file_name} to obtain charge energy information")
            cc_data = None
        return cc_data

    def get_atom_spins(self, cc_data, file):
        atom_spins = getattr(cc_data, "atomspins", None) if cc_data is not None else None
        if atom_spins:
            return atom_spins

        mulliken_spins = self.parse_gaussian_mulliken_spins(file)
        if mulliken_spins is None:
            return {}
        return {"mulliken": mulliken_spins}

    @staticmethod
    def parse_gaussian_mulliken_spins(file):
        spin_tables = []
        with open(file, encoding="utf-8", errors="replace") as handle:
            lines = handle.readlines()

        i = 0
        while i < len(lines):
            if lines[i].strip() != "Mulliken charges and spin densities:":
                i += 1
                continue

            table = []
            i += 1
            while i < len(lines):
                parts = lines[i].split()
                if len(parts) == 4 and parts[0].isdigit():
                    try:
                        table.append(float(parts[3]))
                    except ValueError:
                        break
                elif table:
                    break
                i += 1

            if table:
                spin_tables.append(table)
            continue

        if not spin_tables:
            return None
        return np.array(spin_tables[-1])

    def get_filename(self, fullname):
        return resolve_data_key(fullname, self.data_dict, module_name="CHARGES", logger=self.args.log)

