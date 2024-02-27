######################################################.
#        This file stores the get_df class            #
######################################################.


import pandas as pd
import time
from dftdescp.argument_parser import command_line_args


class get_df:
    """
    Class to create a dataframe of parameters.
    """

    def __init__(self, data_dicts, data_type, substructure=False):
        self.dd = data_dicts
        if data_type == "molecular":
            mol_df = self.get_mol_df()
            return mol_df

        # self.df = False
        # for dict in data_dicts:
        #     dict_df = False
        #     for file_name in dict.keys():
        #         value_column = []
        #         property_column = []
        #         for property in dict[file_name].keys():
        #             value = dict[file_name][property]
        #             value_column.append(value)
        #             property_column.append(property)
        #         if dict_df == False:
        #             data = {"Property": property_column, file_name: value_column}
        #             dict_df = pd.DataFrame(data)
        #         else:
        #             dict_df[file_name] = value_column
        #     if self.df == False:
        #         self.df = dict_df
        #     else:
        #         pd.concat(self.df, dict_df)
        # self.df.to_csv(f'dftdescp_{time.strftime("%Y%m%d-%H%M%S")}_out.csv')
        # return self.df

    # create a df of molecular properties
    def get_mol_df(self):
        mol_df = False
        mol_list = ["opt", "sp_ieea", "ad_ieea"]
        calced_list = self.dd.keys()
        # go through each of the three options for dictionaries of molecular properties
        for category in mol_list:
            # looks to see if these calcs were done
            if category in calced_list:
                dict_df = False
                dict = self.dd[category]
                for file_name in dict.keys():
                    value_column = []
                    property_column = []
                    for property in dict[file_name].keys():
                        value = dict[file_name][property]
                        value_column.append(value)
                        property_column.append(property)
                    if dict_df == False:
                        data = {"Property": property_column, file_name: value_column}
                        dict_df = pd.DataFrame(data)
                    else:
                        dict_df[file_name] = value_column
                if mol_df == False:
                    mol_df = dict_df
                else:
                    pd.concat(mol_df, dict_df)
        return mol_df
