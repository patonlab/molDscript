#!/usr/bin/env python

###############################################.
#          __main__ file for the code         #
###############################################.

from __future__ import absolute_import

import sys
from moldscript import moldscript

# If we are running from a wheel, add the wheel to sys.path
# This allows the usage python pip-*.whl/pip install pip-*.whl

if __package__ != 'moldscript':
    print('moldscript is not installed! Use: pip install moldscript (anywhere, using a terminal) or python setup.py install (from the downloaded /moldscript/moldscript folder).')

if __name__ == '__main__':
    moldscript.main()
    sys.exit()
