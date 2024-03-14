from dftdescp.files import files
from dftdescp.fukui import fukui
from dftdescp.ie_ea import ie_ea
from dftdescp.opt import opt
from dftdescp.nmr import nmr
from dftdescp.nbo import nbo
from dftdescp.substructure import substructure
from dftdescp.get_df import get_df
from dftdescp.argument_parser import command_line_args, dftdescp_version, dftdescp_ref, time_run
import subprocess, sys


header = """
   ▓█████▄   ██████  ▄████▄   ██▀███   ██▓ ██▓███  ▄▄▄█████▓
   ▒██▀ ██▌▒██    ▒ ▒██▀ ▀█  ▓██ ▒ ██▒▓██▒▓██░  ██▒▓  ██▒ ▓▒
   ░██   █▌░ ▓██▄   ▒▓█    ▄ ▓██ ░▄█ ▒▒██▒▓██░ ██▓▒▒ ▓██░ ▒░
   ░▓█▄   ▌  ▒   ██▒▒▓▓▄ ▄██▒▒██▀▀█▄  ░██░▒██▄█▓▒ ▒░ ▓██▓ ░
   ░▒████▓ ▒██████▒▒▒ ▓███▀ ░░██▓ ▒██▒░██░▒██▒ ░  ░  ▒██▒ ░ 
    ▒▒▓  ▒ ▒ ▒▓▒ ▒ ░░ ░▒ ▒  ░░ ▒▓ ░▒▓░░▓  ▒▓▒░ ░  ░  ▒ ░░   
    ░ ▒  ▒ ░ ░▒  ░ ░  ░  ▒     ░▒ ░ ▒░ ▒ ░░▒ ░         ░    
    ░ ░  ░ ░  ░  ░  ░          ░░   ░  ▒ ░░░         ░      
      ░          ░  ░ ░         ░      ░                    
    ░               ░    Paton Research Group, Colorado 2024                              
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
    print("   DFTDESCP v {} {} \n   Citation: {}\n".format(dftdescp_version, time_run, dftdescp_ref))
    print(f'\nArguments passed to program: \n{sys.argv[1:]}\n')
    if args.link:
        # ALL DATA
        all_read = files(calc="link", path=args.path_link)
        if args.opt:
            opt_data = opt(all_read.file_data)
            data_dicts["opt"] = opt_data
        if args.nmr:
            nmr_data = nmr(all_read.file_data)
            data_dicts["nmr"] = nmr_data
        if args.nbo:
            nbo_data = nbo(all_read.file_data)
            data_dicts["nbo"] = nbo_data
        # if args.fukui :  fukui_data = fukui(all_read.file_data) ## I dont think we can access this using link
        # if args.sp_ie_ea:  sp_ie_ea_data = ie_ea(all_read.file_data) ## I dont think we can access this using link
        # if args.ad_ie_ea : ad_ie_ea_data =  ie_ea(all_read.file_data) ## I dont think we can access this using link

    else:
        # OPT
        if args.opt:
            opt_read = files(calc="opt", path=args.path_opt)
            opt_data = opt(opt_read.file_data)
            # print(opt_data.file_data.keys())
            data_dicts["opt"] = opt_data

        # NMR
        if args.nmr:
            nmr_read = files(calc="nmr", path=args.path_nmr)
            nmr_data = nmr(nmr_read.file_data)
            # print(nmr_data.file_data.keys())
            data_dicts["nmr"] = nmr_data
        # NBO
        if args.nbo:
            nbo_read = files(calc="nbo", path=args.path_nbo)
            nbo_data = nbo(nbo_read.file_data)
            # print(nbo_data.file_data.keys())
            data_dicts["nbo"] = nbo_data

        # FUKUI
        if args.fukui:
            fukui_read = files(calc="fukui", path=args.path_fukui)
            fukui_data = fukui(fukui_read.file_data)
            # print(fukui_data.file_data.keys())
            data_dicts["fukui"] = fukui_data

        # SP IE & EA
        if args.sp_ie_ea:
            sp_ie_ea_read = files(calc="sp_ie_ea", path=args.path_sp_ie_ea)
            sp_ie_ea_data = ie_ea(sp_ie_ea_read.file_data)
            # print(sp_ie_ea_data.file_data.keys())
            data_dicts["sp_ieea"] = sp_ie_ea_data

        # AD IE & EA
        if args.ad_ie_ea:
            ad_ie_ea_read = files(calc="ad_ie_ea", path=args.path_ad_ie_ea)
            ad_ie_ea_data = ie_ea(ad_ie_ea_read.file_data)
            # print(ad_ie_ea_data.file_data.keys())
            data_dicts["ad_ieea"] = ad_ie_ea_data
    
    if args.substructure != "":
        substructure_read = files(calc="substructure", path=args.path_opt)
        substructure_data = substructure(substructure_read.file_data, args.substructure)
        #print(
        #    substructure_data.file_data[
        #        "/Users/shreesowndarya/github/dftdecsp/tests/QCALC/success/Ac4_rdkit_conf_1.log"
         #   ][args.substructure]["index"]
        #)
        atom_df = get_df(data_dicts, 'atom', substructure= substructure_data.file_data)
        if args.nbo: bond_df = get_df(nbo_data, 'bond', substructure= substructure_data.file_data)
    
    else:
        atom_df = get_df(data_dicts, 'atom')
        if args.nbo: bond_df = get_df(nbo_data, 'bond')
    mol_df = get_df(data_dicts, 'molecular')
   
    




    # # Creates the parameterizer class and writes a .csv to look at for atom mapping
    # param = parameterizer(args["struc"], args["sdf"], args["txt"])

    # ######Need to do something about this step for sure########
    # atom_labels = {"log_name": "log_name", 0: "C4", 1: "C1", 2: "O3", 3: "H5", 4: "O2"}
    # ###############################################

    # # generate the atom mapped df and writes it to a .csv
    # param.generate_df(atom_labels)

    # # gets the parameters and writes them to a .csv
    # param.get_params(args["skip_list"])


if __name__ == "__main__":
    checks()
    main()
