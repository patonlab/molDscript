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
        #print(self.dd)
        if data_type == "molecular":
            mol_df = self.get_mol_df()
            self.mol_df = mol_df

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
        mol_df = pd.DataFrame()
        mol_list = ["opt", "sp_ieea", "ad_ieea"]
        calced_list = list(self.dd.keys())
        # go through each of the three options for dictionaries of molecular properties
        for category in mol_list:
            # looks to see if these calcs were done
            if category in calced_list:
                dict = self.dd[category].file_data

                start = False
                if category == 'opt':
                    for file_name in dict.keys():
                        for title in dict[file_name].keys():

                            basename = self.file_base(file_name)
                            final_dict = dict[file_name][title]

                            properties = list(final_dict.keys())

                            properties.insert(0, 'File')
                            if start == False:
                                dict_df = {k: [] for k in properties}
                                start = True
                            for property in properties:
                                if property == 'File':
                                    dict_df[property].append(basename)
                                else:
                                    dict_df[property].append(final_dict[property])
                elif category == 'sp_ieea':
                     start = False
                     for file_name in dict.keys():
                            basename = self.file_base(file_name)
                            final_dict = dict[file_name]
                            ie = final_dict['ie']['E']
                            ea = final_dict['ea']['E']
                            if start == False:
                                dict_df = {k: [] for k in ['File', 'sp_ie','sp_ea']}
                                start=True
                            dict_df['File'].append(basename)
                            dict_df['sp_ie'].append(ie)
                            dict_df['sp_ea'].append(ea)
                elif category == 'ad_ieea':
                     start = False
                     for file_name in dict.keys():
                            basename = self.file_base(file_name)
                            final_dict = dict[file_name]
                            ie = final_dict['ie']['E']
                            ea = final_dict['ea']['E']
                            if start == False:
                                dict_df = {k: [] for k in ['File', 'ad_ie','ad_ea']}
                                start=True
                            dict_df['File'].append(basename)
                            dict_df['ad_ie'].append(ie)
                            dict_df['ad_ea'].append(ea)
                dict_df = pd.DataFrame(dict_df)
                if mol_df.empty:
                    mol_df = dict_df
                else:
                    mol_df = mol_df.merge(dict_df,how='left', on='File')
        mol_df.to_csv('mol_df.csv')
        return mol_df
    def file_base(self, string):
        try:
            int(string[-1])
        except:
            pass
        else:
            return string
        for i in string[::-1]:
            try:
                int(i)
            except:
                pass
            else:
                lastidx = string.rfind(i) + 1

                break
        startidx = string.rfind('/') +1

        return string[startidx:lastidx]