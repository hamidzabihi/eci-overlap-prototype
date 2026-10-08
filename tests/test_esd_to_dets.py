"""
Unit tests for the ESD -> detstring mapper.

These tests use synthetic FragmentOccupation objects to validate
the mapping logic, the electron bookkeeping, and the dets file
round-trip. They do NOT require SHARC, wfoverlap, or any QM code.
"""

from pathlib import Path
import numpy as np
import pytest

from src.esd_to_dets import (
    FragmentOccupation,
    esd_to_detstring,
    fragment_to_occupation,
    total_charge,
    total_electrons,
    total_nmo,
    write_dets_file,
    read_dets_file,
    detstring_to_symbol_counts,
)


# ------------------------------------------------------------ fragment unit

def test_empty_fragment():
    frag = FragmentOccupation(label="A", Z=0, S=0, N=1, n_mo=4)
    assert fragment_to_occupation(frag) == ["e", "e", "e", "e"]


def test_single_alpha():
    frag = FragmentOccupation(label="A", Z=0, S=1, N=1, n_mo=4, alpha_occ=[1])
    assert fragment_to_occupation(frag) == ["e", "a", "e", "e"]


def test_single_beta():
    frag = FragmentOccupation(label="A", Z=0, S=1, N=1, n_mo=4, beta_occ=[2])
    assert fragment_to_occupation(frag) == ["e", "e", "b", "e"]


def test_doubly_occupied():
    frag = FragmentOccupation(label="A", Z=0, S=0, N=1, n_mo=4,
                              alpha_occ=[0], beta_occ=[0])
    assert fragment_to_occupation(frag) == ["d", "e", "e", "e"]


def test_mixed_occupation():
    frag = FragmentOccupation(label="A", Z=0, S=0, N=1, n_mo=5,
                              alpha_occ=[0, 2], beta_occ=[0, 3])
    assert fragment_to_occupation(frag) == ["d", "e", "a", "b", "e"]


def test_out_of_range_occupation_raises():
    with pytest.raises(ValueError):
        FragmentOccupation(label="A", Z=0, S=0, N=1, n_mo=3, alpha_occ=[5])


# ------------------------------------------------------------ ESD mapping

def test_esd_two_fragments():
    """Two BODIPY-like fragments, both neutral singlet ground states."""
    bd = FragmentOccupation(label="BD", Z=0, S=0, N=1, n_mo=4,
                            alpha_occ=[0, 1], beta_occ=[0, 1])
    ba = FragmentOccupation(label="BA", Z=0, S=0, N=1, n_mo=3,
                            alpha_occ=[0], beta_occ=[0])
    det = esd_to_detstring([bd, ba])
    assert det == "dd" + "ee" + "d" + "ee"
    assert len(det) == 7
    counts = detstring_to_symbol_counts(det)
    assert counts == {"d": 3, "a": 0, "b": 0, "e": 4}


def test_esd_charge_transfer():
    """BD+ (doublet) and BA- (doublet) — charge transfer ESD."""
    bd = FragmentOccupation(label="BD", Z=1, S=1, N=1, n_mo=3,
                            alpha_occ=[0], beta_occ=[])
    ba = FragmentOccupation(label="BA", Z=-1, S=1, N=1, n_mo=3,
                            alpha_occ=[], beta_occ=[0])
    det = esd_to_detstring([bd, ba])
    assert det == "ae" + "e" + "be" + "e"
    counts = detstring_to_symbol_counts(det)
    assert counts == {"d": 0, "a": 1, "b": 1, "e": 4}


def test_total_electrons_and_charge():
    bd = FragmentOccupation(label="BD", Z=0, S=0, N=1, n_mo=2,
                            alpha_occ=[0], beta_occ=[0])
    ba = FragmentOccupation(label="BA", Z=-1, S=1, N=1, n_mo=2,
                            alpha_occ=[], beta_occ=[0])
    frags = [bd, ba]
    assert total_electrons(frags) == (1, 2)
    assert total_charge(frags) == -1
    assert total_nmo(frags) == 4


# ------------------------------------------------------------ dets round-trip

def test_dets_roundtrip(tmp_path):
    """Write a dets file and read it back; coefficients must match."""
    detstrings = ["ddee", "aebe", "eebe"]
    coeffs = np.array([
        [0.979083342437, -0.094807515471, 0.001],
        [-0.122637656388, -0.663224542162, 0.01],
    ])  # shape (nstate=2, nsd=3)
    out = tmp_path / "dets"
    write_dets_file(out, detstrings, coeffs, nmo=4)

    nstate, nmo, dets_back, coeffs_back = read_dets_file(out)
    assert nstate == 2
    assert nmo == 4
    assert dets_back == detstrings
    np.testing.assert_allclose(coeffs_back, coeffs, atol=1e-10)


def test_dets_file_header(tmp_path):
    """The first line must be 'nstate nMO nSD'."""
    out = tmp_path / "dets"
    write_dets_file(out, ["dd"], np.array([[1.0]]), nmo=2)
    header = out.read_text().splitlines()[0]
    assert header == "1 2 1"


def test_dets_mismatched_shapes_raise(tmp_path):
    with pytest.raises(ValueError):
        write_dets_file(
            tmp_path / "dets",
            ["dd", "ee"],
            np.array([[1.0]]),  # only 1 det, but 2 detstrings
            nmo=2,
        )


def test_dets_wrong_detstring_length_raises(tmp_path):
    with pytest.raises(ValueError):
        write_dets_file(
            tmp_path / "dets",
            ["dd"],           # length 2
            np.array([[1.0]]),
            nmo=4,            # but says 4
        )
