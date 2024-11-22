#####################################################.
#      This file contains the argument parser        #
#####################################################.

import os, time, getopt, sys
from pathlib import Path
from moldscript.utils import format_lists, Logger

moldscript_version = "0.1"
time_run = time.strftime("%Y/%m/%d %H:%M:%S", time.localtime())
moldscript_ref = "King, J.; Sowndarya, S. V. S.; Paton, R. S. MOLDSCRIPT: A General-Purpose Workflow for Quantum Chemical Molecular Descriptors"

# default variables
var_dict = {
    "struc": "",
    "program": "gaussian",
    "varfile": None,
    "command_line": False,
    "verbose": True,
    "opt": False,
    "volume": False,
    "vall": False,
    "spc": False,
    "nmr": False,
    "nbo": False,
    "charges": False,
    "fmo": False,
    "fukui_neutral": False,
    "fukui_oxidized": False,
    "fukui_reduced": False,
    "ad_oxidized": False,
    "ad_reduced": False,
    "skip_list": [],
    "link": False,
    "boltz": False,
    "min_max": False,
    "temp":298.15,
    "cut":0.95,
    "syllables": 1,
    "substructure": "",
    "varfile" : '',
    "spc_program": 'gaussian',
    "radius": 3,
    "output": ""
}


# part for using the options in a script or jupyter notebook
class options_add:
    pass

def load_arguments_from_file(filename):
    """
    Load arguments from a .txt file where arguments are in the format 'arg : value'.
    """
    kwargs = {}
    with open(filename, 'r') as file:
        for line in file:
            # Skip empty lines and comments
            if not line.strip() or line.startswith("#"):
                continue
            
            # Split the line into argument and value
            arg, value = map(str.strip, line.split(":", 1))
            
            # Convert the value to the appropriate type
            if value.lower() == "true":
                value = True
            elif value.lower() == "false":
                value = False
            elif value.isdigit():
                value = int(value)
            else:
                try:
                    value = float(value)
                except ValueError:
                    pass  # Keep value as a string if it can't be converted
            
            kwargs[arg] = value
    return kwargs

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
        "boltz",
        "min_max",
        "vall",
        "volume"
    ]
    list_args = ["skip_list"]
    int_args = ["syllables"]
    float_args = [
        "temp",
        "cut"
    ]
    str_args = [
        "output",
        "struct",
        "program",
        "varfile",
        "opt",
        "spc",
        "nmr",
        "nbo",
        "charges",
        "fukui_neutral",
        "fukui_oxidized",
        "fukui_reduced",
        "ad_reduced",
        "ad_oxidized",
        "link",
        "substructure",
        "varfile",
        "radius",
        "fmo"
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
                f"o  MOLDSCRIPT v {moldscript_version} is installed correctly! For more information about the available options, see the documentation in XXX"
            )
            sys.exit()
        else:
            # this converts the string parameters to lists
            if arg_name == "varfile" and value:
                # Load arguments from the .txt file
                file_kwargs = load_arguments_from_file(value)
                kwargs.update(file_kwargs)
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

    # Check if 'opt' keyword is provided
    if 'opt' not in kwargs:
        print("Error: The 'opt' keyword is required.")
        sys.exit(1)

    # Second, load all the default variables as an "add_option" object
    args = load_variables(kwargs, "command")

    if len(args.skip_list) != 0:
        for prop in args.skip_list:
            vars(args)[prop] = False

    return args


def load_variables(kwargs, moldscript_module, create_dat=True):
    """
    Load default and user-defined variables
    """
    # first, load default values and options manually added to the function
    self = set_options(kwargs)

    if moldscript_module != "command":
        self.initial_dir = Path(os.getcwd())

        error_setup = False

        # start a log file to track the  module
        if create_dat:
            if moldscript_module == "FILES":
                logger_1 = "FILES"

            if moldscript_module == "OPT":
                logger_1 = "OPT"
            
            elif moldscript_module == "SPC":
                logger_1 = "SPC"

            elif moldscript_module == "NBO":
                logger_1 = "NBO"

            elif moldscript_module == "NMR":
                logger_1 = "NMR"

            elif moldscript_module == "IE_EA":
                logger_1 = "IE_EA"

            elif moldscript_module == "FUKUI":
                logger_1 = "FUKUI"

            if moldscript_module == "SUBSTRUCTURE":
                logger_1 = "SUBSTRUCTURE"
            if moldscript_module == "FMO":
                logger_1 = "FMO"
            if not error_setup:
                if not self.command_line:
                    self.log = Logger(
                        f"{self.initial_dir}/MOLDSCRIPT", logger_1, verbose=self.verbose
                    )
                else:
                    # prevents errors when using command lines and running to remote directories
                    path_command = Path(f"{os.getcwd()}")
                    self.log = Logger(
                        f"{path_command}/MOLDSCRIPT", logger_1, verbose=self.verbose
                    )

                self.log.write_only(
                    f"   MOLDSCRIPT v {moldscript_version} {time_run} \n   Citation: {moldscript_ref}\n"
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
                        f"Command line used in MOLDSCRIPT: python -m moldscript {cmd_print}\n"
                    )

            if error_setup:
                # this is added to avoid path problems in jupyter notebooks
                self.log.finalize()
                os.chdir(self.initial_dir)
                sys.exit()
    return self
