from dftdescp.parameterizer import parameterizer
from dftdescp.files import files
from dftdescp.argument_parser import command_line_args
import subprocess, sys

def checks():
    # this is a dummy import just to warn the user if Open babel is not installed
    try:
        command_run_1 = ["obabel", "-H"]
        subprocess.run(command_run_1, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except FileNotFoundError:
        print("x  Open Babel is not installed! You can install the program with 'conda install -c conda-forge openbabel'")
        sys.exit()
    try: 
        from rdkit.Chem import AllChem as Chem
    except ModuleNotFoundError:
        print("x  RDKit is not installed! You can install the program with 'conda install -c conda-forge rdkit'")
        sys.exit()

def main():
    # This chunk parses the CLI arguments and load user-defined arguments from command line
    args = command_line_args()  

    if args.link:
        # ALL DATA
        all_read = files(calc='link',path=args.path_link)
        print(all_read.file_data.keys())
        # nmr_data = nmr(all_read.file_data)
        # nbo_data = nmr(all_read.file_data)
        # fukui_data = nmr(all_read.file_data)
        # sp_ie_ea_data = nmr(all_read.file_data)
        # ad_ie_ea_data = nmr(all_read.file_data)

    else:
        # NMR
        if args.nmr:
            nmr_read = files(calc='nmr',path=args.path_nmr)
            print(nmr_read.file_data.keys())
            # nmr_data = nmr(nmr_read.file_data)

        # NBO
        if args.nbo:
            nbo_read = files(calc='nbo',path=args.path_nbo)
            print(nbo_read.file_data.keys())
            # nbo_data = nmr(nbo_read.file_data)

        # FUKUI
        if args.fukui:
            fukui_read = files(calc='fukui',path=args.path_fukui)
            print(fukui_read.file_data.keys())
            # fukui_data = nmr(fukui_read.file_data)

        # SP IE & EA
        if args.sp_ie_ea:
            sp_ie_ea_read = files(calc='sp_ie_ea',path=args.path_sp_ie_ea)
            print(sp_ie_ea_read.file_data.keys())
            # sp_ie_ea_data = nmr(sp_ie_ea_read.file_data)

        # AD IE & EA
        if args.ad_ie_ea:
            ad_ie_ea_read = files(calc='ad_ie_ea',path=args.path_ad_ie_ea)
            print(ad_ie_ea_read.file_data.keys())
            # ad_ie_ea_data = nmr(sp_ie_ea_read.file_data)
    
    

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