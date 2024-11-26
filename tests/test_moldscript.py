#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from conftest import datapath
from moldscript.files import files
from moldscript.opt import opt


@pytest.mark.parametrize("path, program, species, energy, smiles", [
    ('arbr/opt', 'gaussian', 'arbr31_wb97xd', -2996.631169953559, '[H]c1c([H])c2c(c([H])c1Br)C([H])([H])C([H])([H])C2=O'),
])

def test_arbr_opt(path, program, species, energy, smiles):
    path = datapath(path)
    data_dicts = {}

    opt_read = files("opt", path, data_dicts, program = program)
    opt_data = opt(opt_read.file_data, data_dicts, program = program)
    data_dicts = opt_data.file_data
    
    precision = 6 # 6 decimal places for energies
    
    # molecular parameters
    assert round(data_dicts[species]['mol']['scfenergy'], precision) == round(energy, precision)
    assert data_dicts[species]['mol']['smiles'] == smiles


    

