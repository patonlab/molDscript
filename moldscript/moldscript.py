import subprocess, sys
from moldscript.files import files
from moldscript.fukui import fukui
from moldscript.opt import opt
from moldscript.spc import spc
from moldscript.nmr import nmr
from moldscript.nbo import nbo
from moldscript.substructure import substructure
from moldscript.get_df import get_df
from moldscript.min_max import min_max
from moldscript.sterics import sterics
from moldscript.charges import charges
from moldscript.lowe import lowe
from moldscript.argument_parser import (
    command_line_args,
    moldscript_version,
    moldscript_ref,
    time_run,
)
from moldscript.boltz import boltz
from moldscript.fmo import fmo
import time

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



def main():
    # This chunk parses the CLI arguments and load user-defined arguments from command line
    args = command_line_args()
    tstart = time.time()

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
            opt_read = files("opt", args.opt, data_dicts, args.suffix_opt)
            opt_data = opt(opt_read.file_data, data_dicts)
            data_dicts = opt_data.file_data
        
        # SPC
        if args.spc:
            spc_read = files(calc="spc", path=args.spc, data_dict=data_dicts, suffix=args.suffix_spc)
            spc_data = spc(spc_read.file_data, data_dicts)
            data_dicts = spc_data.file_data
        
        # Charges
        if args.charges:
            chg_read = files(calc="charges", path=args.charges, data_dict=data_dicts, suffix=args.suffix_charges)
            chg_data = charges(chg_read.file_data, data_dicts)
            data_dicts = chg_data.file_data
        
        # FMO
        if args.fmo:
            fmo_read = files(calc="fmo", path=args.fmo, data_dict=data_dicts, suffix=args.suffix_fmo)
            fmo_data = fmo(fmo_read.file_data, data_dicts)
            data_dicts = fmo_data.file_data
        
        # NMR
        if args.nmr:
            nmr_read = files("nmr", args.nmr, data_dicts, args.suffix_nmr)
            nmr_data = nmr(nmr_read.file_data, data_dicts)
            data_dicts = nmr_data.file_data

        # NBO
        if args.nbo:
            nbo_read = files("nbo", args.nbo, data_dicts, args.suffix_nbo)
            nbo_data = nbo(nbo_read.file_data, data_dicts)
            data_dicts = nbo_data.file_data

        # FUKUI
        if args.fukui_neutral and args.fukui_reduced and args.fukui_oxidized:
            print('FUKUI PATH', [args.fukui_neutral, args.fukui_reduced, args.fukui_oxidized])
            fukui_read = files(calc="fukui", data_dict=data_dicts, path=[args.fukui_neutral, args.fukui_reduced, args.fukui_oxidized], suffix= [args.suffix_fukui_neutral, args.suffix_fukui_reduced, args.suffix_fukui_oxidized])
            fukui_data = fukui(fukui_read.file_data, data_dicts)
            data_dicts = fukui_data.data_dict
    
    if args.substructure != "":
        substructure_read = files(data_dict=data_dicts, calc="substructure", path=args.opt)
        data_dicts = substructure(substructure_read.file_data, data_dicts, args.substructure).file_data
    
    if args.volume != False or args.vall != False:
        data_dicts = sterics(opt_read.file_data, data_dicts, args.volume, args.vall, args.radius).dd
            
    df_getter = get_df(data_dicts, substructure=args.substructure, prefix = args.output, bond_filter=args.no_bond_filter)
    
    if args.boltz:
        boltz(temp=args.temp, prefix=args.output, energies = df_getter.energies)
    
    if args.min_max:
        min_max(temp=args.temp, cut=args.cut,  prefix=args.output, energies = df_getter.energies)
    if args.lowe:
        lowe(prefix=args.output, energies = df_getter.energies)
    
    tfin = time.time()
    print(F"\n\tMolDscript finished running in {round(tfin - tstart, 2)} seconds")
if __name__ == "__main__":
    checks()
    main()