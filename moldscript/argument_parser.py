#####################################################.
#      This file contains the argument parser        #
#####################################################.

import os, time, getopt, sys, shlex
from moldscript.utils import format_lists, Logger, build_log_path

moldscript_version = "0.1"
time_run = time.strftime("%Y/%m/%d %H:%M:%S", time.localtime())
moldscript_ref = "King, J.; Sowndarya, S. V. S.; Paton, R. S. MOLDSCRIPT: A General-Purpose Workflow for Quantum Chemical Molecular Descriptors"

# default variables
var_dict = {
    "struc": "",
    "varfile": None,
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
    "link": False,
    "boltz": False,
    "min_max": False,
    "lowe" : False,
    "temp":298.15,
    "cut":0.95,
    "syllables": 1,
    "substructure": "",
    "varfile" : '',
    "radius": 3,
    "output": "",
    "ml" : True,
    "no_bond_filter": False,
    "suffix_opt": "",
    "suffix_spc": "",
    "suffix_nbo": "",
    "suffix_nmr": "",
    "suffix_fukui_neutral": "",
    "suffix_fukui_reduced": "",
    "suffix_fukui_oxidized": "",
    "suffix_charges": "",
    "suffix_fmo": "",
    "no_mol" : False,
    'no_atom' : False,
    'no_bond' : False,
    'mol_vector' : False,

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
        "boltz",
        "min_max",
        "lowe",
        "vall",
        "volume",
        'ml',
        'no_bond_filter',
        'no_mol',
        'no_atom',
        'no_bond',
        'mol_vector'
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
        "varfile",
        "opt",
        "spc",
        "nmr",
        "nbo",
        "charges",
        "fukui_neutral",
        "fukui_oxidized",
        "fukui_reduced",
        "link",
        "substructure",
        "varfile",
        "radius",
        "fmo",
        "opt_suffix",
        "spc_suffix",
        "nbo_suffix",
        "nmr_suffix",
        "fukui_neutral_suffix",
        "fukui_reduced_suffix",
        "fukui_oxidized_suffix",
        "charges_suffix",
        "fmo_suffix"
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
    # Second, load all the default variables as an "add_option" object
    args = load_variables(kwargs, "command")

    return args


def load_variables(kwargs, moldscript_module, create_dat=True):
    """
    Load default and user-defined variables
    """
    # first, load default values and options manually added to the function
    self = set_options(kwargs)

    self.log = Logger.silent()

    if moldscript_module != "command" and create_dat:
        module_codes = {
            "FILES": "FILES",
            "OPT": "OPT",
            "SPC": "SPC",
            "NBO": "NBO",
            "NMR": "NMR",
            "IE_EA": "IE_EA",
            "FUKUI": "FUKUI",
            "SUBSTRUCTURE": "SUBSTRUCTURE",
            "FMO": "FMO",
        }
        logger_code = module_codes.get(moldscript_module, moldscript_module)
        log_path = build_log_path(self.output, logger_code)
        self.log = Logger(log_path, verbose=self.verbose)
        self.log.write_only(
            f"   MOLDSCRIPT v {moldscript_version} {time_run} \n   Citation: {moldscript_ref}\n"
        )
        command_line = shlex.join(["python", "-m", "moldscript", *sys.argv[1:]])
        self.log.write(f"Command line used in MOLDSCRIPT: {command_line}")

    return self


