"""ORCA parser regression tests.

These guard against drift in cclib's ORCA parser and in the per-module
scrapers (`orca_nmr_shielding`, `bondorders`, `npa_data`) when cclib is
upgraded. The numerical expectations come from running each module on
the fixtures at the time of writing; if cclib changes how it surfaces
ORCA output, these tests will fail and force a deliberate update.

Fixtures: tests/orca/{QCALC,NMR,NBO} — arbr_01..28 r2scan3c conformers.
"""
import pytest
from conftest import fixturepath
from moldscript.files import files
from moldscript.opt import opt
from moldscript.nmr import nmr
from moldscript.nbo import nbo


ORCA_OPT = fixturepath("orca/QCALC")
ORCA_NMR = fixturepath("orca/NMR")
ORCA_NBO = fixturepath("orca/NBO")
SPECIES = "arbr_01"


def _run_opt():
    data_dicts = {}
    opt_read = files("opt", ORCA_OPT, data_dicts, "r2scan3c")
    return opt(opt_read.file_data, data_dicts, create_dat=False).file_data


def test_orca_opt_smiles_and_energy():
    dd = _run_opt()
    assert SPECIES in dd
    assert dd[SPECIES]["mol"]["smiles"] == "[H]c1c([H])c2c(c([H])c1Br)C([H])([H])OC2=O"
    # SCF energy comes from ORCA via cclib (eV → hartree).
    # 4 decimal places ≈ ±0.0001 Ha tolerance is well below physical
    # meaningfulness and well above any parser-precision drift.
    assert round(dd[SPECIES]["mol"]["scfenergy"], 4) == round(-3032.3852, 4)


def test_orca_nmr_shielding_tensors():
    dd = _run_opt()
    nmr_read = files("nmr", ORCA_NMR, dd, "nmr")
    dd = nmr(nmr_read.file_data, dd, create_dat=False).file_data

    shifts = dd[SPECIES]["atom"]["nmr_shielding"]
    expected = [
        -58.854, 11.328, 116.409, 112.032, 28.198, 53.951, 29.943,
        1990.968, 46.355, 51.937, 51.48, 26.797, 26.798, 24.273, 24.158, 23.896,
    ]
    assert len(shifts) == len(expected)
    for i, (actual, want) in enumerate(zip(shifts, expected)):
        assert round(actual, 2) == round(want, 2), f"shift {i}: got {actual}, want {want}"


def test_orca_nbo_charges_and_bond_orders():
    dd = _run_opt()
    # ORCA NBO files use a different suffix (`nbo`) on top of the
    # `r2scan3c_wB97MV` token — `get_filename` in utils.py walks back
    # one `_<token>` at a time to recover the `arbr_01` opt key.
    nbo_read = files("nbo", ORCA_NBO, dd, "nbo")
    dd = nbo(nbo_read.file_data, dd, create_dat=False).file_data

    charges = dd[SPECIES]["atom"]["natural_charge"]
    bond_orders = dd[SPECIES]["atom"]["bond_orders"]
    expected_charges = [
        -0.54664, 0.7705, -0.46272, -0.09114, -0.00024, -0.23819, -0.05381, 0.07235,
        -0.23927, -0.13755, -0.18312, 0.20259, 0.20261, 0.23259, 0.23315, 0.23889,
    ]
    expected_bonds = [
        2.0764, 3.8505, 2.1915, 3.8322, 4.0055, 3.9548, 4.0225, 1.2033,
        3.958, 3.9509, 3.9976, 0.9609, 0.9609, 0.9481, 0.9474, 0.9444,
    ]
    for i, (actual, want) in enumerate(zip(charges, expected_charges)):
        assert round(actual, 3) == round(want, 3), f"charge {i}: got {actual}, want {want}"
    for i, (actual, want) in enumerate(zip(bond_orders, expected_bonds)):
        assert round(actual, 2) == round(want, 2), f"bond_order {i}: got {actual}, want {want}"
    # Sanity check: ORCA natural charges should sum to total charge (0 here).
    assert round(sum(charges), 2) == 0.0


@pytest.mark.skip(
    reason="xTB outputs are not recognized by cclib's ccopen() (only by ccread "
    "via openbabel fallback, which returns an empty ccData). The hand-rolled "
    "xTB CPU-time branch in opt.py never runs because initiate_data_dict() "
    "uses the strict cclib path. xTB opt is therefore broken end-to-end."
)
def test_xtb_opt_smoke():
    pass
