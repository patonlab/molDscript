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
from moldscript.argument_parser import command_line_args, moldscript_version, moldscript_ref, time_run
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
    print("   MOLDSCRIPT v {} {} \n   {}\n".format(moldscript_version, time_run, moldscript_ref))
    print(f'   Arguments passed to program: \n   {sys.argv[1:]}\n')

    if args.path_opt is not None and not args.opt: args.opt = True
    if args.path_nmr is not None and not args.nmr: args.nmr = True
    if args.path_nbo is not None and not args.nbo: args.nbo = True
    if args.path_fukui is not None and not args.fukui: args.fukui = True
    if args.path_sp_ie_ea is not None and not args.sp_ie_ea: args.sp_ie_ea = True

    if args.link:
        # ALL DATA
        all_read = files(calc="link", path=args.path_link, program=args.program)
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
            
            opt_read = files(calc="opt", path=args.path_opt, 
                             program=args.program, suffix_opt=args.suffix_opt)
            opt_data = opt(opt_read.file_data, data_dicts,
                           program=args.program, volume=args.volume)
            data_dicts = opt_data.file_data
            
        # SPC
        if args.spc:
            spc_read = files(calc="spc", path=args.path_spc,
                            suffix_spc=args.suffix_spc,
                            program=args.spc_program)
            spc_data = spc(spc_read.file_data, data_dicts,
                           spc_program=args.spc_program)
            data_dicts = spc_data.file_data
            
        # NMR
        if args.nmr:
            nmr_read = files(calc="nmr", path=args.path_nmr, 
                             suffix_nmr=args.suffix_nmr,
                             program=args.program)
            nmr_data = nmr(nmr_read.file_data, data_dicts,
                           program=args.program)
            data_dicts = nmr_data.file_data

        # NBO
        if args.nbo:
            nbo_read = files(calc="nbo", path=args.path_nbo,
                            suffix_nbo=args.suffix_nbo,
                            program=args.program)
            nbo_data = nbo(nbo_read.file_data, data_dicts,
                           program=args.program)
            data_dicts["nbo"] = nbo_data

        # FUKUI
        if args.fukui:
            fukui_read = files(calc="fukui", path=args.path_fukui, 
                               suffix_fukui=args.suffix_fukui, 
                               suffix_fred=args.suffix_fred,
                               suffix_fox=args.suffix_fox,
                               program=args.program)
            fukui_data = fukui(fukui_read.file_data,
                               program=args.program)
            data_dicts["fukui"] = fukui_data

        # SP IE & EA
        if args.sp_ie_ea:
            print(args.path_sp_ie_ea, args.suffix_sp_ie,args.suffix_sp_ea)

            sp_ie_ea_read = files(calc="sp_ie_ea", path=args.path_sp_ie_ea, 
                                  suffix_sp_ie=args.suffix_sp_ie,
                                  suffix_sp_ea=args.suffix_sp_ea,
                                  program=args.program)
            sp_ie_ea_data = ie_ea(sp_ie_ea_read.file_data,
                                  program=args.program)
            data_dicts["sp_ieea"] = sp_ie_ea_data

        # AD IE & EA
        if args.ad_ie_ea:
            ad_ie_ea_read = files(calc="ad_ie_ea", path=args.path_ad_ie_ea,
                                  suffix_ad_ie=args.suffix_ad_ie,
                                  suffix_ad_ea=args.suffix_ad_ea,
                                  program=args.program)
            ad_ie_ea_data = ie_ea(ad_ie_ea_read.file_data,
                                  program=args.program)
            data_dicts["ad_ieea"] = ad_ie_ea_data
    
    if args.substructure != "":
        substructure_read = files(calc="substructure", path=args.path_opt)
        substructure_data = substructure(substructure_read.file_data, args.substructure)

        if args.fukui or args.nmr or args.nbo or args.volume: atom_df = get_df(data_dicts, 'atom', substructure= substructure_data.file_data, program=args.program)
        if args.nbo or args.opt : bond_df = get_df(data_dicts, 'bond', substructure= substructure_data.file_data, nbo_suffix=args.suffix_nbo, program=args.program)
    
    else:
        if args.fukui or args.nmr or args.nbo or args.volume: atom_df = get_df(data_dicts, 'atom', program=args.program, volume=args.volume)
        if args.nbo or args.opt: bond_df = get_df(data_dicts, 'bond', program=args.program)

    if args.opt: 
        mol_df = get_df(data_dicts, 'molecular', nbo_suffix=args.suffix_nbo, program=args.program)
    if args.boltz: 
        boltz(temp=args.temp, spc=args.spc, syllables=args.syllables)
    if args.min_max: 
        min_max(temp=args.temp, cut=args.cut, spc=args.spc, syllables=args.syllables)

if __name__ == "__main__":
    checks()
    main()
