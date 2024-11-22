######################################################.
#          This file stores functions used           #
#                in multiple modules                 #
######################################################.

import os
import ast
from pathlib import Path
import glob

k_B_hartree = 3.1668114e-6  # hartree/K
J_TO_AU = 4.184 * 627.509541 * 1000.0  # UNIT CONVERSION
eV_to_hartree = 0.0367493

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

def add_cpu_times(file_data):
    ''' add cpu times for all files'''
    total_cpu = datetime.timedelta(0)
    for filename in file_data.keys():
        total_cpu += file_data[filename]['cpu_time']
    
    return total_cpu

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
    if program == 'xtb':
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
def get_filename(fullname, dd):
    flist = list(dd.keys())
    tempname = fullname
    for i in range(fullname.count("_")):
        try:
            findex = flist.index(tempname)
            keyname = flist[findex]
            return keyname
        except:
            tempname = tempname.rsplit("_", 1)[0]
            print(tempname)
    print(
        f"Error processing file {fullname}. Ensure consistent naming as described in the docs."
    )
    raise SystemExit


