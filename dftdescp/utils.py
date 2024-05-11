######################################################.
#          This file stores functions used           #
#                in multiple modules                 #
######################################################.

import os
import ast
from pathlib import Path
import glob
import yaml

# class for logging
class Logger:
    """
    Class that wraps a file object to abstract the logging.
    """

    # Class Logger to write output to a file
    def __init__(self, filein, append, suffix="dat", verbose=True):
        if verbose:
            self.log = open(f"{filein}_{append}.{suffix}", "w")
        else:
            self.log = ''

    def write(self, message):
        """
        Appends a newline character to the message and writes it into the file.

        Parameters
        ----------
        message : str
           Text to be written in the log file.
        """
        try:
            self.log.write(f"{message}\n")
        except AttributeError:
            pass
        print(f"{message}")

    def write_only(self, message):
        """
        Appends a newline character to the message and writes it into the file, but does not print it.

        Parameters
        ----------
        message : str
           Text to be written in the log file.
        """
        try:
            self.log.write(f"{message}\n")
        except AttributeError:
            pass

    def finalize(self):
        """
        Closes the file
        """
        try:
            self.log.close()
        except AttributeError:
            pass

def format_lists(value):
    '''
    Transforms strings into a list
    '''

    if not isinstance(value, list):
        try:
            value = ast.literal_eval(value)
        except (SyntaxError, ValueError):
            # this line fixes issues when using "[X]" or ["X"] instead of "['X']" when using lists
            value = value.replace('[',']').replace(',',']').replace("'",']').split(']')
            while('' in value):
                value.remove('')
    return value

# load paramters from yaml file
def load_from_yaml(self):
    """
    Loads the parameters for the calculation from a yaml if specified. Otherwise
    does nothing.
    """

    txt_yaml = f"\no  Importing DFTDESCP parameters from {self.varfile}"
    error_yaml = False
    # Variables will be updated from YAML file
    try:
        if os.path.exists(self.varfile):
            if os.path.basename(Path(self.varfile)).split('.')[1] in ["yaml", "yml", "txt"]:
                with open(self.varfile, "r") as file:
                    try:
                        param_list = yaml.load(file, Loader=yaml.SafeLoader)
                    except yaml.scanner.ScannerError:
                        txt_yaml = f'\nx  Error while reading {self.varfile}. Edit the yaml file and try again (i.e. use ":" instead of "=" to specify variables)'
                        error_yaml = True
        if not error_yaml:
            for param in param_list:
                if hasattr(self, param):
                    if getattr(self, param) != param_list[param]:
                        setattr(self, param, param_list[param])

    except UnboundLocalError:
        txt_yaml = f"\nx  The specified yaml file containing parameters {self.varfile} was not found or the extension is not compatible ('.yaml', '.yml' or '.txt')! Also, make sure that the params file is in the folder where you are running the code."

    return self, txt_yaml

def get_files(value, program):
    if program == 'gaussian':
        if value[-1]=='/':
            value = value[:-1]
        if (
        Path(f"{value}").exists()
        and os.getcwd() not in f"{value}"
        ):
            list_of_val = glob.glob(f"{os.getcwd()}/{value}/*.log")
        else:
            list_of_val = glob.glob(value)
        return list_of_val
    if program == 'orca':
        if value[-1]=='/':
            value = value[:-1]
        if (
        Path(f"{value}").exists()
        and os.getcwd() not in f"{value}"
        ):
            list_of_val = glob.glob(f"{os.getcwd()}/{value}/*.out")

        else:
            list_of_val = glob.glob(value)

        return list_of_val
def find_nth(haystack: str, needle: str, n: int) -> int:
    start = haystack.find(needle)
    while start >= 0 and n > 1:
        start = haystack.find(needle, start+len(needle))
        n -= 1
    return start



