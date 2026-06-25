#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pathlib import Path

import pandas as pd
import pytest
from conftest import datapath
from moldscript.files import files
from moldscript.MLIP import mlip
from moldscript.get_df import get_df
from moldscript.opt import opt
from moldscript.nmr import nmr
from moldscript.nbo import nbo
from moldscript.fmo import fmo
from moldscript.charges import charges
from moldscript.fukui import fukui
from moldscript.utils import eV_to_hartree

@pytest.mark.parametrize("opt_path, species,  smiles, suffix",[
    ('arbr/opt',  'arbr31_wb97xd',  '[H]c1c([H])c2c(c([H])c1Br)C([H])([H])C([H])([H])C2=O', '',)
])

def test_opt(opt_path,  species,  smiles, suffix):
    path = datapath(opt_path)
    data_dicts = {}

    opt_read = files("opt", path, data_dicts, suffix)
    opt_data = opt(opt_read.file_data, data_dicts)
    data_dicts = opt_data.file_data
    
    # molecular parameters
    precision = 6 # 6 decimal places for energies
    assert data_dicts[species]['mol']['smiles'] == smiles


@pytest.mark.parametrize("opt_path, nmr_path, species, shifts, nmr_suffix", [
    ('arbr/opt', 'arbr/nmr',  'arbr31_wb97xd', [-257.5005, -20.2833, 149.8784, 159.2040, 24.3340, 
     52.8746, 34.6443, 2046.8050, 52.1611, 57.5798, 45.7539, 29.2925, 29.2924, 28.7646, 28.7645, 24.0512, 24.1362, 23.8194], 'nmr'),
])

def test_nmr(opt_path, nmr_path,  species, shifts, nmr_suffix):
    path = datapath(opt_path)
    data_dicts = {}
    path = datapath(nmr_path)
    nmr_read = files("nmr", path, data_dicts, nmr_suffix)
    nmr_data = nmr(nmr_read.file_data, data_dicts, nmr_suffix)
    data_dicts = nmr_data.file_data
    
    # molecular parameters
    precision = 2 # 2 decimal places for chemical shifts                      
    for i, shift in enumerate(shifts):
        print(data_dicts.keys())
        assert round(data_dicts[species]['atom']['nmr_shielding'][i], precision) == round(shift, precision)


@pytest.mark.parametrize("opt_path, nbo_path,  species, natural_charges, bond_orders, nbo_suffix", [
    ('arbr/opt', 'arbr/popn',  'arbr31_wb97xd', [-0.52472, 0.57494, -0.50777, -0.42749, 
    0.05447, -0.2471, -0.0508, 0.06753, -0.24629, -0.13454, -0.18304, 0.2398, 0.2398, 0.22423, 0.22423, 0.2296, 0.23156, 0.2356],
    [2.0669, 3.8883, 3.9071, 3.916, 4.0067, 3.9575, 4.0203, 1.2074, 3.9584, 3.9514, 4.0012, 0.9433, 0.9433, 0.9509, 0.9509, 0.949, 0.9477, 0.9458], 'popn'),
])

def test_nbo(opt_path, nbo_path,  species, natural_charges, bond_orders, nbo_suffix):
    path = datapath(opt_path)
    data_dicts = {}

    path = datapath(nbo_path)
    nbo_read = files("nbo", path, data_dicts, nbo_suffix)
    nbo_data = nbo(nbo_read.file_data, data_dicts, )
    data_dicts = nbo_data.file_data
    
    # molecular parameters
    precision = 2 # 2 decimal places for charges and bond orders                      
    for i, charge in enumerate(natural_charges):
        assert round(data_dicts[species]['atom']['natural_charge'][i], precision) == round(charge, precision)
    assert round(sum(data_dicts[species]['atom']['natural_charge']), precision) == round(0, precision)
    for i, order in enumerate(bond_orders):
        assert round(data_dicts[species]['atom']['bond_orders'][i], precision) == round(order, precision)


@pytest.mark.parametrize("opt_path, species, dipole, homo, lumo, quadrupole_moment_trace, fmo_suffix", [
    ('arbr/opt', 'arbr31_wb97xd', 2.536, -0.33291, -0.00473, -234.10379, ''),
])

def test_fmo(opt_path, species, dipole, homo, lumo, quadrupole_moment_trace, fmo_suffix):
    path = datapath(opt_path)
    data_dicts = {}
    fmo_read = files("fmo", path, data_dicts, fmo_suffix)
    fmo_data = fmo(fmo_read.file_data, data_dicts, fmo_suffix)
    data_dicts = fmo_data.file_data
    
    # molecular parameters
    precision = 2 # 2 decimal places for dipole        
    assert round(data_dicts[species]['mol']['dipole'], precision) == round(dipole, precision)
    assert round(data_dicts[species]['mol']['quadrupole_moment_trace'], precision) == round(quadrupole_moment_trace, precision)
    
    precision = 6 # 6 decimal places for energies        
    assert round(data_dicts[species]['mol']['HOMO'] * eV_to_hartree, precision) == round(homo, precision)
    assert round(data_dicts[species]['mol']['LUMO'] * eV_to_hartree, precision) == round(lumo, precision)

    # test operations: HOMO-LUMO_gap, chemical_potential, global_electrophilicity, global_nucleophilicity
    hl_gap = lumo - homo
    hl_sum = homo + lumo
    chemical_potential = hl_sum / 2
    global_electrophilicity = chemical_potential ** 2 / hl_gap / 2
    global_nucleophilicity = 1 / global_electrophilicity
    
    precision = 4 # 4 decimal places    
    assert round(data_dicts[species]['mol']['HOMO-LUMO_gap'] * eV_to_hartree, precision) == round(hl_gap, precision)
    assert round(data_dicts[species]['mol']['chemical_potential'] * eV_to_hartree, precision) == round(chemical_potential, precision)
    assert round(data_dicts[species]['mol']['global_electrophilicity'] * eV_to_hartree, precision) == round(global_electrophilicity, precision)
    assert round(data_dicts[species]['mol']['global_nucleophilicity'] / eV_to_hartree, precision) == round(global_nucleophilicity, precision)


@pytest.mark.parametrize("opt_path, species, apt_charges, charges_suffix", [
    ('arbr/opt',  'arbr31_wb97xd', [-0.813669,  1.084436, -0.094902,  0.06463 ,  0.119114, -0.205275,
        0.450111, -0.26189 , -0.244025,  0.105473, -0.406945,  0.018934,
        0.018939, -0.012211, -0.012213,  0.061321,  0.065093,  0.063078], ''),
])

def test_charges(opt_path,  species, apt_charges, charges_suffix):
    path = datapath(opt_path)
    data_dicts = {}
    chg_read = files("charges", path, data_dicts, charges_suffix)
    chg_data = charges(chg_read.file_data, data_dicts, charges_suffix)
    data_dicts = chg_data.file_data

    # molecular parameters
    precision = 2 # 2 decimal places for charges and bond orders                      
    for i, charge in enumerate(apt_charges):
        assert round(data_dicts[species]['atom']['apt_charge'][i], precision) == round(charge, precision)


@pytest.mark.parametrize("opt_path, fukui_neutral_path, fukui_oxidized_path, fukui_reduced_path,  fukui_neutral_suffix, fukui_oxidized_suffix, fukui_reduced_suffix, species, oxidized_charges, reduced_charges", [
    ('arbr/opt', 'arbr/popn', 'arbr/fukui_ox', 'arbr/fukui_red', 'popn', 'ox', 'red', 'arbr31_wb97xd', [-0.07481, 0.55092, -0.42669, -0.43292, 0.09900, -0.24454, 
    -0.00150, 0.17203, -0.22431, -0.11532, -0.16606, 0.28381, 0.28382, 0.26198, 0.26197, 0.25935, 0.26189, 0.25138], [-0.67645, 0.44630,-0.48153, -0.41843, 0.03051, -0.31967, -0.17854, -0.05180, -0.23203, -0.26912, -0.22398, 0.19031, 0.19031, 0.19610, 0.19610, 0.19775, 0.19804, 0.20612]),
])

def test_fukui(opt_path, fukui_neutral_path, fukui_oxidized_path, fukui_reduced_path, fukui_neutral_suffix, fukui_oxidized_suffix, fukui_reduced_suffix, species, oxidized_charges, reduced_charges):
    path = datapath(opt_path)
    data_dicts = {}
    fukui_read = files(
                calc="fukui",
                data_dict=data_dicts,
                path=[datapath(fukui_neutral_path), datapath(fukui_reduced_path), datapath(fukui_oxidized_path)],
                suffix=[fukui_neutral_suffix, fukui_reduced_suffix, fukui_oxidized_suffix])
    fukui_data = fukui(fukui_read.file_data, data_dicts)
    data_dicts = fukui_data.data_dict

    # molecular parameters - here checking the charges and fukui indices sum to the right values
    precision = 2 # 2 decimal places for charges and bond orders                           
    for i, charge in enumerate(oxidized_charges):
        assert round(data_dicts[species]['atom']['oxidized_natural_charges'][i], precision) == round(charge, precision)
    
    for i, charge in enumerate(reduced_charges):
        assert round(data_dicts[species]['atom']['reduced_natural_charges'][i], precision) == round(charge, precision)
                
    assert round(sum(data_dicts[species]['atom']['oxidized_natural_charges']), precision) == round(1, precision)
    assert round(sum(data_dicts[species]['atom']['reduced_natural_charges']), precision) == round(-1, precision)
    assert round(sum(data_dicts[species]['atom']['fminus']), precision) == round(1, precision)
    assert round(sum(data_dicts[species]['atom']['fplus']), precision) == round(1, precision)
    assert round(sum(data_dicts[species]['atom']['frad']), precision) == round(1, precision)  


def build_arbr141_parser():
    return mlip(
        neutral=datapath("arbr/mlip/neutral"),
        reduced=datapath("arbr/mlip/reduced"),
        oxidized=datapath("arbr/mlip/oxidized"),
    )


def test_mlip_parses_arbr141_dictionary_shape():
    parser = build_arbr141_parser()

    assert "CPU_time" in parser.file_data
    assert "arbr141_wb97xd" in parser.file_data

    entry = parser.file_data["arbr141_wb97xd"]
    assert set(entry.keys()) == {"mol", "atom", "bond", "CPU_time"}
    assert entry["bond"]["bond_length"].shape == (15, 15)
    assert len(entry["atom"]["atomnos"]) == 15


def test_mlip_arbr141_numeric_values():
    parser = build_arbr141_parser()
    entry = parser.file_data["arbr141_wb97xd"]

    assert entry["mol"]["scfenergy"] == pytest.approx(-3176.0668656171874, rel=0, abs=1e-9)
    assert entry["mol"]["vertical_ie"] == pytest.approx(0.33734708984375, rel=0, abs=1e-12)
    assert entry["mol"]["vertical_ea"] == pytest.approx(0.04105585859375, rel=0, abs=1e-12)
    assert entry["mol"]["dipole_x"] == pytest.approx(-1.1098721027374268, rel=0, abs=1e-12)
    assert entry["mol"]["dipole_y"] == pytest.approx(-0.11217078566551208, rel=0, abs=1e-12)
    assert entry["mol"]["dipole_z"] == pytest.approx(-0.17991302907466888, rel=0, abs=1e-12)
    assert entry["mol"]["dipole"] == pytest.approx(1.1299411788329683, rel=0, abs=1e-12)

    assert entry["atom"]["charges_neutral"].sum() == pytest.approx(0.0, abs=1e-6)
    assert entry["atom"]["fplus"].sum() == pytest.approx(1.0, abs=1e-6)
    assert entry["atom"]["fminus"].sum() == pytest.approx(1.0, abs=1e-6)
    assert entry["atom"]["frad"].sum() == pytest.approx(1.0, abs=1e-6)


def test_mlip_get_df_writes_csvs(tmp_path):
    parser = build_arbr141_parser()

    prefix = str(tmp_path / "mlip_")
    get_df(parser.file_data, prefix=prefix)

    mol_csv = Path(f"{prefix}molecule_level.csv")
    atom_csv = Path(f"{prefix}atom_level.csv")
    bond_csv = Path(f"{prefix}bond_level.csv")

    assert mol_csv.exists()
    assert atom_csv.exists()
    assert bond_csv.exists()

    mol_df = pd.read_csv(mol_csv)
    atom_df = pd.read_csv(atom_csv)
    bond_df = pd.read_csv(bond_csv)

    assert list(mol_df["filename"]) == ["arbr141_wb97xd"]
    assert {"filename", "atom_index", "atom_type", "fplus", "fminus", "frad"}.issubset(atom_df.columns)
    assert {"filename", "atom1_idx", "atom2_idx", "bond_length"}.issubset(bond_df.columns)