from parameterizer import parameterizer
import pandas as pd


class nmr(parameterizer):
    def __init__(self, a_list):
        super().__init__(self)  # Will need to match these inputs
        self.atom_list = a_list

    def get_nmr(
        dataframe, a_list
    ):  # a function to get the nbo for all atoms (a_list, form ["C1", "C4", "O2"]) in a dataframe that contains file name and atom number
        nmr_df = pd.DataFrame(columns=[])  # define an empty df to place results in

        for index, row in dataframe.iterrows():  # iterate over the dataframe

            # if True:
            try:  # try to get the data
                atom_list = []
                for new_a in a_list:
                    new_atom = row[
                        str(new_a)
                    ]  # the atom number (i.e. 16) to add to the list is the df entry of this row for the labeled atom (i.e.) "C1")
                    atom_list.append(
                        str(new_atom)
                    )  # append that to atom_list to make a list of the form [16, 17, 29]
                log_file = row["log_name"]  # read file name from df
                filecont, error = get_filecont(
                    log_file
                )  # read the contents of the log file
                if error != "":
                    print(error)
                    row_i = {}
                    for a in range(0, len(a_list)):
                        entry = {"NMR_shift_" + str(a_list[a]): "no data"}
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
                    for a in range(0, len(a_list)):
                        entry = {"NMR_shift_" + str(a_list[a]): "no data"}
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
                nmrout = get_specdata(atom_list, nmr)  # revisit
                # print(nmrout)

                # this adds the data from the nboout into the new property df
                row_i = {}
                for a in range(0, len(a_list)):
                    entry = {"NMR_shift_" + str(a_list[a]): nmrout[a]}
                    row_i.update(entry)
                nmr_dataframe = nmr_dataframe.append(row_i, ignore_index=True)
            except:
                print("****Unable to acquire NMR shifts for:", row["log_name"], ".log")
                row_i = {}
                for a in range(0, len(a_list)):
                    entry = {"NMR_shift_" + str(a_list[a]): "no data"}
                    row_i.update(entry)
                nmr_dataframe = nmr_dataframe.append(row_i, ignore_index=True)
        print("NMR function has completed for", a_list)
        return pd.concat([dataframe, nmr_dataframe], axis=1)
