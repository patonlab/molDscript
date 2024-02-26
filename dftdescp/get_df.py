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

    def __init__(self, data_dicts):
        self.df = False
        for dict in data_dicts:
            dict_df = False
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
            if self.df == False:
                self.df = dict_df
            else:
                pd.concat(self.df, dict_df)
        self.df.to_csv(f'dftdescp_{time.strftime("%Y%m%d-%H%M%S")}_out.csv')
        return self.df
