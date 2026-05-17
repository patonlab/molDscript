#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import pytest

try:
    import moldscript
    BASEPATH = os.path.join(moldscript.__path__[0])
except ImportError:
    here = os.path.dirname(os.path.abspath(__file__))
    BASEPATH = os.path.normpath(os.path.join(here, '..', 'moldscript'))

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))


def datapath(path):
    """Path under moldscript/examples/ — Gaussian fixtures shipped with the package."""
    return os.path.join(BASEPATH, 'examples', path)


def fixturepath(path):
    """Path under tests/ — fixtures kept outside the installed package
    (ORCA, xTB, etc.). Renamed away from `testpath` to avoid pytest's
    collection picking it up as a test function."""
    return os.path.join(TESTS_DIR, path)

