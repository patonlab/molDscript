#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from conftest import datapath
from moldscript.files import files
from moldscript.opt import opt


@pytest.mark.parametrize("path, program, energy", [
    ('arbr/opt/', 'gaussian', 298.15),
])

def test_arbr_opt(path, program, energy):
    path = datapath(path)
    data_dicts = {}

    opt_read = files("opt", path, data_dicts, program = program)
    #print("printing opt_read", opt_read.file_data)
    opt_data = opt(opt_read.file_data, data_dicts, program = program)
    #data_dicts = opt_data.file_data
    
    precision = 6 # if temp == 298.15 else 4e-4
    
    print(data_dicts)
    #assert e == round(bbe.scf_energy, precision)
    

