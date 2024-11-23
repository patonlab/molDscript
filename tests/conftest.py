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


def datapath(path):
    return os.path.join(BASEPATH, 'tests', path)

