#!/usr/bin/env python

###############################################.
#          __main__ file for the code         #
###############################################.

from __future__ import absolute_import

import sys
from dftdescp import moldscript

# If we are running from a wheel, add the wheel to sys.path
# This allows the usage python pip-*.whl/pip install pip-*.whl

if __package__ != 'dftdescp':
    print('dftdescp is not installed! Use: pip install dftdescp (anywhere, using a terminal) or python setup.py install (from the downloaded /dftdescp/dftdescp folder).')

if __name__ == '__main__':
    moldscript.main()
    sys.exit()
