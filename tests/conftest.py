#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import pytest

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))


def datapath(path):
    """Path under tests/data/ — Gaussian fixtures used by the test suite.

    These used to live at moldscript/examples/ but were moved out of the
    installed package to avoid shipping ~2 MB of .log files to every
    user who pip-installs moldscript.
    """
    return os.path.join(TESTS_DIR, 'data', path)


def fixturepath(path):
    """Path under tests/ — fixtures kept outside the installed package
    (ORCA, xTB, etc.). Renamed away from `testpath` to avoid pytest's
    collection picking it up as a test function."""
    return os.path.join(TESTS_DIR, path)

