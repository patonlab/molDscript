from rdkit import Chem
import pandas as pd
import get_properties_functions as gp


class parameterizer:
    """
    This class parameterizes substructures of molecules in a dataset for ML applications based on DFT optputs.

    Attributes
    ----------
    prelim_df: pandas df
        dataframe before atom labeling
    mapped_df: pandas df
        dataframe after atom labeling that gets parameters added to it

    Methods
    -------
    def generate_df(self, atom_labels)
        generates and saves the atom mapped df
    def get_params(self, skip_list)
        Calculates the parameters not in the skip list and adds them to the mapped_df attribute then saves this df
    """

    def __init__(self, struc, sdf, txt):
        """
        Creates an instance of the parameterizer class.

        Parameters
        ----------
        struc: str
            The smarts substructure you are looking at
        sdf: path
            Path to open the .sdf file
        txt: path
            Path to the .txt file



        """

        substructure = Chem.MolFromSmarts(struc)
        # generate a list of molecules using RDkit
        all_compounds = Chem.SDMolSupplier(sdf, removeHs=False)
        # molecules.sdf is generated with the instructions above
        # it is a single sdf that contains the structures/atom numbers etc. for every molecule you will analyze
        # uses RDKit to search for the substructure in each compound you will analyze
        atoms = []
        for molecule in all_compounds:
            if molecule is not None:
                submatch = molecule.GetSubstructMatches(
                    substructure
                )  # find substructure
                matchlist = list(
                    [item for sublist in submatch for item in sublist]
                )  # list of zero-indexed atom numbers
                match_idx = [
                    x + 1 for x in matchlist
                ]  # this line changes from 0-indexed to 1-indexed (for Gaussian)
                atoms.append(
                    match_idx
                )  # append 1-indexed list to atoms (a list of lists)

        # this loop extracts log names from log_ids and splits them to the desired format
        filenames = open(txt, "r")  # generate this with instruction above
        # it is a text file that contains the file name for every molecule you will analyze
        list_of_filenames = [
            (line.strip()).split() for line in filenames
        ]  # list of the file names (each of which includes all conformers)
        list_of_files = []
        for filename in list_of_filenames:
            file = filename[0].split(".")
            list_of_files.append(file[0])
        filenames.close()

        # put the atom numbers for the substructure for each log file into a dataframe
        prelim_df = pd.DataFrame(atoms)
        prelim_df.insert(0, column="log_name", value=list_of_files)
        self.prelim_df = prelim_df

    def get_filecont(log):  # gets the entire job output
        error = ""  # default unless "normal termination" is in file
        an_error = True
        with open(log + ".log") as f:
            loglines = f.readlines()
        for line in loglines[::-1]:
            if "Normal termination" in line:
                an_error = False
            if an_error:
                error = "****Failed or incomplete jobs for " + log + ".log"
        return (loglines, error)

    def generate_df(self, atom_labels):
        """
        Generate the atom mapped df

        Parameters
        ----------
        atom_labels: dict
            prelimary df column labels based on user input from gview"""

        self.mapped_df = self.prelim_df.rename(columns=atom_labels)
        self.mapped_df.to_csv("mapped_df.csv")
        return self.mapped_df

    def get_params(self, skip_list):
        ##For now ths still uses their get_properties_functions file but I want to transfer those to be functions of this class
        """
        Expand the df with the parameters you are aquiring.

        Paramters
        ---------
        skip_list: list
            List of the methods you do not want to run ex: ['goodvibes_e', 'frontierorbs']
        """

        # ---------------GoodVibes Engergies---------------
        # uses the GoodVibes 2021 Branch (Jupyter Notebook Compatible)
        # calculates the quasi harmonic corrected G(T) and single point corrected G(T) as well as other thermodynamic properties
        # inputs: dataframe, temperature
        if "goodvibes_e" not in skip_list:
            self.df = gp.get_goodvibes_e(self.mapped_df, 298.15)

        # ---------------Frontier Orbitals-----------------
        # E(HOMO), E(LUMO), mu(chemical potential or negative of molecular electronegativity), eta(hardness/softness), omega(electrophilicity index)
        if "frontierorbs" not in skip_list:
            self.mapped_df = gp.get_frontierorbs(self.mapped_df)

        # ---------------Polarizability--------------------
        # Exact polarizability
        if "polarizability" not in skip_list:
            self.mapped_df = gp.get_polarizability(self.mapped_df)

        # There will be way more of these, just three for a placeholder now
        self.mapped_df.to_csv("final_parameterset.csv")
        return self.mapped_df
