#####################################################.
#      This file contains the argument parser        #
#####################################################.

import os, time, getopt, sys
from pathlib import Path
from dftdescp.utils import format_lists, load_from_yaml, Logger

dftdescp_version = "0.1"
time_run = time.strftime("%Y/%m/%d %H:%M:%S", time.localtime())
dftdescp_ref = "XXX"

# default variables
var_dict = {
    "struc": "",
    "program": "gaussian",
    "varfile": None,
    "command_line": False,
    "verbose": True,
    "all": False,
    "opt": False,
    "spc": False,
    "nmr": False,
    "nbo": False,
    "fukui": False,
    "ad_ie_ea": False,
    "sp_ie_ea": False,
    "skip_list": [],
    "link": False,
    "boltz": False,
    "min_max": False,
    "temp":298.15,
    "cut":0.95,
    "syllables": 1,
    "substructure": "",
    
    "path_opt": None,
    "path_spc": None,
    "path_nmr": None,
    "path_nbo": None,
    "path_fukui": None,
    "path_ad_ie_ea": None,
    "path_sp_ie_ea": None,
    "path_link": None,
    
    "suffix_opt": None,
    "suffix_spc": "spc",
    "suffix_nmr": "SP_NMR",
    "suffix_nbo": "SP_NBO",
    "suffix_fukui": "SP_neutral",
    "suffix_fred": "SP_reduced",
    "suffix_fox": "SP_oxidized",
    "suffix_ad_ie": "AD_IE",
    "suffix_ad_ea": "AD_EA",
    "suffix_sp_ie": "SP_IE",
    "suffix_sp_ea": "SP_EA",
    "spc_program": 'gaussian',
}


# part for using the options in a script or jupyter notebook
class options_add:
    pass


def set_options(kwargs):
    # set default options and options provided
    options = options_add()
    # dictionary containing default values for options

    for key in var_dict:
        vars(options)[key] = var_dict[key]
    for key in kwargs:
        if key in var_dict:
            vars(options)[key] = kwargs[key]
        elif key.lower() in var_dict:
            vars(options)[key.lower()] = kwargs[key.lower()]
        else:
            print(
                "Warning! Option: [",
                key,
                ":",
                kwargs[key],
                "] provided but no option exists, try the online documentation to see available options for each module.",
            )
    return options


def command_line_args():
    """
    Load default and user-defined arguments specified through command lines. Arrguments are loaded as a dictionary
    """

    # First, create dictionary with user-defined arguments
    kwargs = {}
    available_args = ["help"]
    bool_args = [
        "command_line",
        "all",
        "opt",
        "spc",
        "nmr",
        "nbo",
        "fukui",
        "ad_ie_ea",
        "sp_ie_ea",
        "link",
        "boltz",
        "min_max"
    ]
    list_args = ["skip_list"]
    int_args = ["syllables"]
    float_args = [
        "temp",
        "cut"
    ]
    str_args = [
        "struct",
        "program",
        "varfile",
        "path_opt",
        "path_spc",
        "path_nmr",
        "path_nbo",
        "path_fukui",
        "path_ad_ie_ea",
        "path_sp_ie_ea",
        "path_link",
        "substructure",
        "suffix_spc",
        "suffix_nmr",
        "suffix_nbo",
        "suffix_fukui",
        "suffix_fred",
        "suffix_fox",
        "suffix_ad_ie",
        "suffix_ad_ea",
        "suffix_sp_ie",
        "suffix_sp_ea",
        "suffix_opt",
        "spc_path"
    ]

    for arg in var_dict:
        if arg in bool_args:
            available_args.append(f"{arg}")
        else:
            available_args.append(f"{arg} =")

    try:
        opts, _ = getopt.getopt(sys.argv[1:], "h", available_args)
    except getopt.GetoptError as err:
        print(err)
        sys.exit()

    for arg, value in opts:
        if arg.find("--") > -1:
            arg_name = arg.split("--")[1].strip()
        elif arg.find("-") > -1:
            arg_name = arg.split("-")[1].strip()

        if arg_name in ("h", "help"):
            print(
                f"o  DFTDESCP v {dftdescp_version} is installed correctly! For more information about the available options, see the documentation in XXX"
            )
            sys.exit()
        else:
            # this converts the string parameters to lists
            if arg_name in bool_args:
                value = True
            elif arg_name.lower() in list_args:
                value = format_lists(value)
            elif arg_name.lower() in int_args:
                value = int(value)
            elif arg_name.lower() in float_args:
                value = float(value)
            elif arg_name.lower() in str_args:
                value = str(value)
            elif value == "None":
                value = None
            elif value == "False":
                value = False
            elif value == "True":
                value = True

            kwargs[arg_name] = value

    # Second, load all the default variables as an "add_option" object
    args = load_variables(kwargs, "command")

    # reassinging properties based on all or specific ones
    if args.all:
        vars(args)["opt"] = True
        vars(args)["nmr"] = True
        vars(args)["nbo"] = True
        vars(args)["fukui"] = True
        vars(args)["sp_ie_ea"] = True
        vars(args)["ad_ie_ea"] = True

    if len(args.skip_list) != 0:
        for prop in args.skip_list:
            vars(args)[prop] = False

    return args


def load_variables(kwargs, dftdescp_module, create_dat=True):
    """
    Load default and user-defined variables
    """

    # first, load default values and options manually added to the function
    self = set_options(kwargs)

    # this part loads variables from yaml files (if varfile is used)
    txt_yaml = ""
    if self.varfile is not None:
        self, txt_yaml = load_from_yaml(self)

    if dftdescp_module != "command":
        self.initial_dir = Path(os.getcwd())

        error_setup = False

        # start a log file to track the  module
        if create_dat:
            if dftdescp_module == "FILES":
                logger_1 = "FILES"

            if dftdescp_module == "OPT":
                logger_1 = "OPT"
            
            elif dftdescp_module == "SPC":
                logger_1 = "SPC"

            elif dftdescp_module == "NBO":
                logger_1 = "NBO"

            elif dftdescp_module == "NMR":
                logger_1 = "NMR"

            elif dftdescp_module == "IE_EA":
                logger_1 = "IE_EA"

            elif dftdescp_module == "FUKUI":
                logger_1 = "FUKUI"

            if dftdescp_module == "SUBSTRUCTURE":
                logger_1 = "SUBSTRUCTURE"

            if txt_yaml not in [
                "",
                f"\no  Importing DFTDESCP parameters from {self.varfile}",
                "\nx  The specified yaml file containing parameters was not found! Make sure that the valid params file is in the folder where you are running the code.\n",
            ]:
                self.log = Logger(
                    f"{self.initial_dir}/DFTDESCP", logger_1, verbose=self.verbose
                )
                self.log.write(txt_yaml)
                error_setup = True

            if not error_setup:
                if not self.command_line:
                    self.log = Logger(
                        f"{self.initial_dir}/DFTDESCP", logger_1, verbose=self.verbose
                    )
                else:
                    # prevents errors when using command lines and running to remote directories
                    path_command = Path(f"{os.getcwd()}")
                    self.log = Logger(
                        f"{path_command}/DFTDESCP", logger_1, verbose=self.verbose
                    )

                self.log.write_only(
                    f"   DFTDESCP v {dftdescp_version} {time_run} \n   Citation: {dftdescp_ref}\n"
                )

                if self.command_line:
                    cmd_print = ""
                    cmd_args = sys.argv[1:]
                    for i, elem in enumerate(cmd_args):
                        if elem[0] in ['"', "'"]:
                            elem = elem[1:]
                        if elem[-1] in ['"', "'"]:
                            elem = elem[:-1]
                        if elem != "-h" and elem.split("--")[-1] not in var_dict:
                            if (
                                cmd_args[i - 1].split("--")[-1] in var_dict
                            ):  # check if the previous word is an arg
                                cmd_print += f'"{elem}'
                            if (
                                i == len(cmd_args) - 1
                                or cmd_args[i + 1].split("--")[-1] in var_dict
                            ):  # check if the next word is an arg, or last word in command
                                cmd_print += f'"'
                        else:
                            cmd_print += f"{elem}"
                        if i != len(cmd_args) - 1:
                            cmd_print += " "

                    self.log.write(
                        f"Command line used in DFTDESCP: python -m dftdescp {cmd_print}\n"
                    )

            if error_setup:
                # this is added to avoid path problems in jupyter notebooks
                self.log.finalize()
                os.chdir(self.initial_dir)
                sys.exit()
    return self
