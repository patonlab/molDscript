from moldscript.files import files
from moldscript.fukui import fukui
from moldscript.ie_ea import ie_ea
from moldscript.opt import opt
from moldscript.spc import spc
from moldscript.nmr import nmr
from moldscript.nbo import nbo
from moldscript.substructure import substructure
from moldscript.get_df import get_df
from moldscript.min_max import min_max
from moldscript.argument_parser import (
    command_line_args,
    moldscript_version,
    moldscript_ref,
    time_run,
)
from moldscript.boltz import boltz
import subprocess, sys


header = """
   • ▌ ▄ ·.       ▄▄▌  ·▄▄▄▄  .▄▄ ·  ▄▄· ▄▄▄  ▪   ▄▄▄·▄▄▄▄▄
   ·██ ▐███▪▪     ██•  ██▪ ██ ▐█ ▀. ▐█ ▌▪▀▄ █·██ ▐█ ▄█•██  
   ▐█ ▌▐▌▐█· ▄█▀▄ ██▪  ▐█· ▐█▌▄▀▀▀█▄██ ▄▄▐▀▀▄ ▐█· ██▀· ▐█.▪
   ██ ██▌▐█▌▐█▌.▐▌▐█▌▐▌██. ██ ▐█▄▪▐█▐███▌▐█•█▌▐█▌▐█▪·• ▐█▌·
   ▀▀  █▪▀▀▀ ▀█▄▀▪.▀▀▀ ▀▀▀▀▀•  ▀▀▀▀ ·▀▀▀ .▀  ▀▀▀▀.▀    ▀▀▀ 
                             Paton Research Group, CO 2024                              
"""


def checks():
    # this is a dummy import just to warn the user if Open babel is not installed
    try:
        command_run_1 = ["obabel", "-H"]
        subprocess.run(
            command_run_1, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        )
    except FileNotFoundError:
        print(
            "x  Open Babel is not installed! You can install the program with 'conda install -c conda-forge openbabel'"
        )
        sys.exit()
    try:
        from rdkit.Chem import AllChem as Chem
    except ModuleNotFoundError:
        print(
            "x  RDKit is not installed! You can install the program with 'conda install -c conda-forge rdkit'"
        )


#         sys.exit()


def main():
    # This chunk parses the CLI arguments and load user-defined arguments from command line
    args = command_line_args()


    data_dicts = {}

    print(header)
    print(
        "   MOLDSCRIPT v {} {} \n   {}\n".format(
            moldscript_version, time_run, moldscript_ref
        )
    )
    print(f"   Arguments passed to program: \n   {sys.argv[1:]}\n")

    try:
        check_redox_vars(args.fukui_neutral, args.fukui_reduced, args.fukui_oxidized, 'fukui')
    except ValueError as e:
        print(e)
        sys.exit()

    if args.link:
        # ALL DATA
        all_read = files(calc="link", path=args.link, program=args.program)
        if args.opt:
            opt_data = opt(all_read.file_data, program=args.program)
            data_dicts["opt"] = opt_data
        if args.nmr:
            nmr_data = nmr(all_read.file_data, program=args.program)
            data_dicts["nmr"] = nmr_data
        if args.nbo:
            nbo_data = nbo(all_read.file_data, program=args.program)
            data_dicts["nbo"] = nbo_data

    else:
        # OPT
        if args.opt:

            opt_read = files("opt", args.opt, data_dicts, program=args.program)
            opt_data = opt(
                opt_read.file_data, data_dicts, program=args.program, volume=args.volume
            )
            data_dicts = opt_data.file_data

        # SPC
        if args.spc:
            spc_read = files(
                calc="spc",
                path=args.spc,
                suffix_spc=args.suffix_spc,
                program=args.spc_program,
            )
            spc_data = spc(spc_read.file_data, data_dicts, spc_program=args.spc_program)
            data_dicts = spc_data.file_data

        # NMR
        if args.nmr:
            nmr_read = files("nmr", args.nmr, data_dicts, program=args.program)
            nmr_data = nmr(nmr_read.file_data, data_dicts, program=args.program)
            data_dicts = nmr_data.file_data

        # NBO
        if args.nbo:
            nbo_read = files("nbo", args.nbo, data_dicts, program=args.program)
            nbo_data = nbo(nbo_read.file_data, data_dicts, program=args.program)
            data_dicts = nbo_data.file_data

        # FUKUI
        if args.fukui_neutral and args.fukui_reduced and args.fukui_oxidized:
            fukui_read = files(
                calc="fukui",
                data_dict=data_dicts,
                path=[args.fukui_neutral, args.fukui_reduced, args.fukui_reduced],
                program=args.program,
            )
            fukui_data = fukui(fukui_read.file_data, data_dicts, program=args.program)
            data_dicts = fukui_data.file_data

        # SP IE & EA
        if args.sp_ie_ea:
            print(args.sp_ie_ea, args.suffix_sp_ie, args.suffix_sp_ea)

            sp_ie_ea_read = files(
                calc="sp_ie_ea",
                path=args.sp_ie_ea,
                suffix_sp_ie=args.suffix_sp_ie,
                suffix_sp_ea=args.suffix_sp_ea,
                program=args.program,
            )
            sp_ie_ea_data = ie_ea(
                "sp_ie_ea", sp_ie_ea_read.file_data, data_dicts, program=args.program
            )
            data_dicts = sp_ie_ea_data.file_data

        # AD IE & EA
        if args.ad_ie_ea:
            ad_ie_ea_read = files(
                calc="ad_ie_ea",
                path=args.ad_ie_ea,
                suffix_ad_ie=args.suffix_ad_ie,
                suffix_ad_ea=args.suffix_ad_ea,
                program=args.program,
            )
            ad_ie_ea_data = ie_ea(
                "ad_ie_ea", ad_ie_ea_read.file_data, data_dicts, program=args.program
            )
            data_dicts = ad_ie_ea_data.file_data

    if args.substructure != "":
        substructure_read = files(data_dicts, calc="substructure", path=args.opt)
        data_dicts = substructure(
            substructure_read.file_data, args.substructure
        ).file_data

    df = get_df(data_dicts, program=args.program, substructure=args.substructure)
    if args.boltz:
        boltz(temp=args.temp, spc=args.spc, syllables=args.syllables)
    if args.min_max:
        min_max(temp=args.temp, cut=args.cut, spc=args.spc, syllables=args.syllables)
def check_redox_vars(var1, var2, var3, calc):
    # Count how many variables are strings
    string_count = sum(isinstance(var, str) for var in [var1, var2, var3])
    # Raise an error if only one or two variables are strings
    if 1 <= string_count <= 2:
        raise ValueError(f"Specify a neutral, reduced, AND oxidized {calc} path, or specify none of them.")


if __name__ == "__main__":
    checks()
    main()
