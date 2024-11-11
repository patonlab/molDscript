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
                path=[args.fukui_neutral, args.fukui_reduced, args.fukui_oxidized],
                program=args.program,
            )
            fukui_data = fukui(fukui_read.file_data, data_dicts, program=args.program)
            data_dicts = fukui_data.file_data

        # AD IE & EA
        if args.ad_reduced and args.ad_oxidized:
            ad_ie_ea_read = files(
                calc="ad_ie_ea",
                data_dict = data_dicts,
                path=[args.ad_reduced, args.ad_oxidized],
                program=args.program,
            )
            ad_ie_ea_data = ie_ea(
                "ad_ie_ea", ad_ie_ea_read.file_data, data_dicts, program=args.program
            )
            data_dicts = ad_ie_ea_data.file_data

    if args.substructure != "":
        substructure_read = files(data_dict=data_dicts, calc="substructure", path=args.opt)
        data_dicts = substructure(substructure_read.file_data, data_dicts, args.substructure).file_data

    df = get_df(data_dicts, program=args.program, substructure=args.substructure)
    if args.boltz:
        boltz(temp=args.temp, spc=args.spc, syllables=args.syllables)
    if args.min_max:
        min_max(temp=args.temp, cut=args.cut, spc=args.spc, syllables=args.syllables)



if __name__ == "__main__":
    checks()
    main()
