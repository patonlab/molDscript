######################################################.
#        This file stores the FMO class               #
######################################################.


import sys, os
import time
import datetime
from moldscript.argument_parser import load_variables
from moldscript.utils import (
    format_timedelta,
    get_filename,
    initiate_data_dict,
    progress_iter,
    record_cpu_time,
    safe_parse,
)
import numpy as np

class fmo:
    """
    Class containing all the functions for the FMO module related to output files
    """

    def __init__(self, data, data_dict, create_dat=True,  **kwargs):

        start_time_overall = time.time()
        # load default and user-specified variables
        self.args = load_variables(kwargs, "FMO", create_dat=create_dat)
        self.data = data
        self.data_dict = data_dict
        self.module_cpu_seconds = 0.0
        if self.data_dict == {}:
            self.data_dict = initiate_data_dict(self.data, logger=self.args.log)
        if len(self.data.keys()) == 0:
            self.args.log.write(f"\nx  Could not find files to obtain information for FMO and moment analysis")
            self.args.log.finalize()
            sys.exit()
        else:
            self.file_data = self.get_data()

        if create_dat:
            elapsed_time = round(time.time() - start_time_overall, 2)
            self.args.log.write(f"-- FMO Collection complete in {elapsed_time} seconds\n")

    def get_data(self):

        self.args.log.write(f"-- FMO Collection starting")
        self.module_cpu_seconds = 0.0

        for idx, file_name in progress_iter(self.data.keys(), self.args.log):
            if self.data[file_name].rsplit('.', 1)[1] == 'log':
                self.fmo_program = 'gaussian'
            elif self.data[file_name].rsplit('.', 1)[1] == 'out':
                self.fmo_program = 'orca'
            fmo_data = safe_parse(self.data[file_name], logger=self.args.log, context="FMO and moment data")
            file_name = get_filename(file_name, self.data_dict, logger=self.args.log)
            if idx == 1 and fmo_data is not None:
                try:
                    self.args.log.write(f"   Functional used: {fmo_data.metadata['functional']}")
                    self.args.log.write(f"   Basis set used: {fmo_data.metadata['basis_set']}")
                except (AttributeError, KeyError):
                    pass

            self.args.log.write_only(f"o  Parsing FMO and Moment Data from {os.path.basename(file_name)}")

            self.data_dict[file_name]["mol"]["dipole"] = np.sqrt(np.sum((fmo_data.moments[0] - fmo_data.moments[1]) ** 2, axis=0))
            self.data_dict[file_name]["mol"]["HOMO"] = fmo_data.moenergies[0][fmo_data.homos[0]]
            self.data_dict[file_name]["mol"]["LUMO"] = fmo_data.moenergies[0][fmo_data.homos[0] + 1]

            softness = self.data_dict[file_name]["mol"]["LUMO"]- self.data_dict[file_name]["mol"]["HOMO"]
            self.data_dict[file_name]["mol"]["HOMO-LUMO_gap"] = (softness)
            chemical_potential = (self.data_dict[file_name]["mol"]["LUMO"] + self.data_dict[file_name]["mol"]["HOMO"]) /2
            self.data_dict[file_name]["mol"]["chemical_potential"] = chemical_potential
            glob_electrophilicity = chemical_potential**2 / (2*softness)
            self.data_dict[file_name]["mol"]["global_electrophilicity"] = glob_electrophilicity
            self.data_dict[file_name]["mol"]["global_nucleophilicity"] = 1/glob_electrophilicity

            try:
                quadrupole_moments = fmo_data.moments[2]
                quadrupole_matrix = np.array([
                    [quadrupole_moments[0], quadrupole_moments[1], quadrupole_moments[2]],
                    [quadrupole_moments[1], quadrupole_moments[3], quadrupole_moments[4]],
                    [quadrupole_moments[2], quadrupole_moments[4], quadrupole_moments[5]],
                ])
                trace = np.trace(quadrupole_matrix)
                self.data_dict[file_name]["mol"]["quadrupole_moment_trace"] = trace
            except (AttributeError, IndexError):
                self.data_dict[file_name]["mol"]["quadrupole_moment_trace"] = None

            cpu_times = fmo_data.metadata.get("cpu_time") if fmo_data and hasattr(fmo_data, "metadata") else None
            record_cpu_time(self.data_dict, file_name, self.data[file_name], cpu_times)

        return self.data_dict


