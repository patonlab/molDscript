from parameterizer import parameterizer
import argparse


# This chunk parses the CLI arguments
parser = argparse.ArgumentParser(description="Take in the values for parameterization")
parser.add_argument("struc", help="Smarts for the substructure of interest")
parser.add_argument("sdf", help="Path to .sdf file input")
parser.add_argument("txt", help="Path to .txt file input")
parser.add_argument("skip_list", help="List of parameters to NOT calculate")
args = vars(parser.parse_args())
print(args["skip_list"])
# Creates the parameterizer class and writes a .csv to look at for atom mapping
param = parameterizer(args["struc"], args["sdf"], args["txt"])

######Need to do something about this step for sure########
atom_labels = {"log_name": "log_name", 0: "C4", 1: "C1", 2: "O3", 3: "H5", 4: "O2"}
###############################################

# generate the atom mapped df and writes it to a .csv
param.generate_df(atom_labels)

# gets the parameters and writes them to a .csv
param.get_params(args["skip_list"])
