#####################################################.
#      This file contains the argument parser        #
#####################################################.

import argparse
import os, time, sys, shlex
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


def build_parser():
    """Build the argparse parser. Defaults come from `var_dict` so the CLI,
    the `--varfile` loader, and programmatic instantiation share one source
    of truth.

    Defaults are intentionally suppressed at the argument level
    (`default=argparse.SUPPRESS`) so we can distinguish "user passed this
    flag" from "user did not" — needed for `--varfile` override semantics.
    """
    parser = argparse.ArgumentParser(
        prog="moldscript",
        description=(
            f"MOLDSCRIPT v{moldscript_version} — converts DFT/QM outputs "
            "into molecule/bond/atom-level descriptor CSVs."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    paths = parser.add_argument_group("calculation paths")
    paths.add_argument("--opt", help="optimization output directory")
    paths.add_argument("--spc", help="single-point output directory (overrides opt SCF energies)")
    paths.add_argument("--nmr", help="NMR output directory")
    paths.add_argument("--nbo", help="NBO output directory")
    paths.add_argument("--charges", help="atomic-charge output directory")
    paths.add_argument("--fmo", help="FMO output directory (HOMO/LUMO + dipole/quadrupole)")
    paths.add_argument("--fukui_neutral", help="neutral charge state output directory")
    paths.add_argument("--fukui_reduced", help="reduced charge state output directory")
    paths.add_argument("--fukui_oxidized", help="oxidized charge state output directory")

    suffixes = parser.add_argument_group(
        "filename suffixes (load-bearing for conformer matching across modules)"
    )
    for calc in ("opt", "spc", "nbo", "nmr", "charges", "fmo",
                 "fukui_neutral", "fukui_reduced", "fukui_oxidized"):
        suffixes.add_argument(f"--suffix_{calc}",
            help=f"trailing _<token> in {calc} filenames used to recover the base key")

    sub = parser.add_argument_group("substructure & sterics")
    sub.add_argument("--substructure", help="SMARTS pattern to filter atom/bond descriptors")
    sub.add_argument("--volume", action="store_true",
        help="compute DBSTEP buried volume on substructure match")
    sub.add_argument("--vall", action="store_true",
        help="compute DBSTEP buried volume for every atom")
    sub.add_argument("--radius",
        help="sphere radius for buried volume: float or list literal '[r1,r2,...]'")

    ensemble = parser.add_argument_group("conformer ensemble aggregation")
    ensemble.add_argument("--boltz", action="store_true",
        help="emit Boltzmann-weighted ensemble CSVs")
    ensemble.add_argument("--min_max", action="store_true",
        help="emit min/max/range CSVs (uses --cut)")
    ensemble.add_argument("--lowe", action="store_true",
        help="emit lowest-energy conformer CSVs")
    ensemble.add_argument("--temp", type=float, help="temperature (K) for Boltzmann weights")
    ensemble.add_argument("--cut", type=float,
        help="population cutoff for min_max filtering (default 0.95)")

    out = parser.add_argument_group("output")
    out.add_argument("--output", help="prefix prepended to every generated filename")
    out.add_argument("--no_bond_filter", action="store_true",
        help="disable the default bond_order>=0.1 / bond_length<=3 filter on bond CSV")
    out.add_argument("--no_mol", action="store_true", help="skip molecule_level.csv")
    out.add_argument("--no_atom", action="store_true", help="skip atom_level.csv")
    out.add_argument("--no_bond", action="store_true", help="skip bond_level.csv")
    out.add_argument("--mol_vector", action="store_true", help="emit molecule-level vector encoding")

    parser.add_argument("--varfile",
        help="path to a `key : value` defaults file; CLI flags override file values")

    # Suppress defaults so the namespace only contains explicitly-passed args.
    # Defaults are applied later via var_dict + varfile + CLI precedence.
    for action in parser._actions:
        if action.dest != "help":
            action.default = argparse.SUPPRESS

    return parser


def command_line_args():
    """Parse CLI args with argparse, layer them on top of `var_dict` defaults
    plus optional `--varfile` overrides, and return a populated `options_add`
    namespace via `load_variables`.

    Precedence (low → high): var_dict defaults → varfile entries → CLI flags.
    """
    parser = build_parser()
    cli_args = vars(parser.parse_args())  # only contains keys actually passed

    kwargs = {}
    varfile_path = cli_args.get("varfile")
    if varfile_path:
        kwargs.update(load_arguments_from_file(varfile_path))
    kwargs.update(cli_args)

    return load_variables(kwargs, "command")


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


