from parameterizer import parameterizer
import pandas as pd
import re


class nmr(parameterizer):
    def __init__(self, a_list):
        super().__init__(self)  # Will need to match these inputs
        self.atom_list = a_list

    def get_nmr(self):

        # these are at the start of the get_properties_functions and are needed for this module
        nmrstart_pattern = " SCF GIAO Magnetic shielding tensor (ppm):\n"
        nmrend_pattern = re.compile("End of Minotr F.D.")
        nmrend_pattern_os = re.compile("g value of the free electron")

        nmr_dataframe = pd.DataFrame(
            columns=[]
        )  # define an empty df to place results in

        for index, row in self.mapped_df.iterrows():  # iterate over the dataframe

            # if True:
            try:  # try to get the data
                atom_list = []
                for new_a in self.a_list:
                    new_atom = row[
                        str(new_a)
                    ]  # the atom number (i.e. 16) to add to the list is the df entry of this row for the labeled atom (i.e.) "C1")
                    atom_list.append(
                        str(new_atom)
                    )  # append that to atom_list to make a list of the form [16, 17, 29]
                log_file = row["log_name"]  # read file name from df
                filecont, error = self.get_filecont(
                    log_file
                )  # read the contents of the log file
                if error != "":
                    print(error)
                    row_i = {}
                    for a in range(0, len(self.a_list)):
                        entry = {"NMR_shift_" + str(self.a_list[a]): "no data"}
                        row_i.update(entry)
                    nmr_dataframe = nmr_dataframe.append(row_i, ignore_index=True)
                    continue

                # determining the locations/values for start and end of NMR section
                start, end, i = 0, 0, 0
                if nmrstart_pattern in filecont:
                    start = filecont.index(nmrstart_pattern) + 1
                    for i in range(start, len(filecont), 1):
                        if nmrend_pattern.search(
                            filecont[i]
                        ) or nmrend_pattern_os.search(filecont[i]):
                            end = i
                            break
                if start == 0:
                    error = (
                        "****no NMR data found in file: "
                        + str(row["log_name"])
                        + ".log"
                    )
                    print(error)
                    row_i = {}
                    for a in range(0, len(self.a_list)):
                        entry = {"NMR_shift_" + str(self.a_list[a]): "no data"}
                        row_i.update(entry)
                    nmr_dataframe = nmr_dataframe.append(row_i, ignore_index=True)
                    continue

                atoms = int(
                    (end - start) / 5
                )  # total number of atoms in molecule (there are 5 lines generated per atom)
                nmr = []
                for atom in range(atoms):
                    element = str.split(filecont[start + 5 * atom])[1]
                    shift_s = str.split(filecont[start + 5 * atom])[4]
                    nmr.append([element, shift_s])
                # atom_list = ["1", "2", "3"]
                nmrout = get_specdata(
                    atom_list, nmr
                )  # Need to figure out what this does - Jake
                # print(nmrout)

                # this adds the data from the nboout into the new property df
                row_i = {}
                for a in range(0, len(self.a_list)):
                    entry = {"NMR_shift_" + str(self.a_list[a]): nmrout[a]}
                    row_i.update(entry)
                nmr_dataframe = nmr_dataframe.append(row_i, ignore_index=True)
            except:
                print("****Unable to acquire NMR shifts for:", row["log_name"], ".log")
                row_i = {}
                for a in range(0, len(self.a_list)):
                    entry = {"NMR_shift_" + str(self.a_list[a]): "no data"}
                    row_i.update(entry)
                nmr_dataframe = nmr_dataframe.append(row_i, ignore_index=True)
        print("NMR function has completed for", self.a_list)
        complete_df = pd.concat([self.mapped_df, nmr_dataframe], axis=1)
        complete_df.to_csv("nmr_parameters.csv")
