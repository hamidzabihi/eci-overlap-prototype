"""
Unit tests for the ESD overlap product formula.
"""

import numpy as np
import pytest

from src.esd_overlap import (
    ESDSpec,
    esd_overlap,
    esd_overlap_matrix,
    fragment_overlaps_from_matrix,
    merge_fragment_overlaps,
)


# ------------------------------------------------------------ ESDSpec

def test_esd_spec_basic():
    esd = ESDSpec(site_states=(("BD", 0), ("BA", 0)))
    assert esd.fragments == ("BD", "BA")
    assert esd.state_of("BD") == 0
    assert esd.state_of("BA") == 0


def test_esd_spec_empty_raises():
    with pytest.raises(ValueError):
        ESDSpec(site_states=())


def test_esd_spec_missing_fragment_raises():
    esd = ESDSpec(site_states=(("BD", 0),))
    with pytest.raises(KeyError):
        esd.state_of("BA")


# ------------------------------------------------------------ esd_overlap

def test_esd_overlap_identity_single_fragment():
    """A single fragment with 1.0 overlap gives ESD overlap 1.0."""
    esd = ESDSpec(site_states=(("BD", 0),))
    O = {("BD", 0, 0): 1.0}
    assert esd_overlap(esd, esd, O) == 1.0


def test_esd_overlap_two_fragments_product():
    """ESD overlap is the product of fragment overlaps."""
    esd_a = ESDSpec(site_states=(("BD", 0), ("BA", 0)))
    esd_b = ESDSpec(site_states=(("BD", 0), ("BA", 0)))
    O = {
        ("BD", 0, 0): 0.9,
        ("BA", 0, 0): 0.8,
    }
    assert abs(esd_overlap(esd_a, esd_b, O) - 0.9 * 0.8) < 1e-15


def test_esd_overlap_mixed_states():
    """Different state indices multiply the right fragment overlaps."""
    esd_a = ESDSpec(site_states=(("BD", 1), ("BA", 0)))
    esd_b = ESDSpec(site_states=(("BD", 0), ("BA", 2)))
    O = {
        ("BD", 1, 0): 0.5,
        ("BA", 0, 2): 0.7,
    }
    assert abs(esd_overlap(esd_a, esd_b, O) - 0.35) < 1e-15


def test_esd_overlap_missing_fragment_key_raises():
    esd_a = ESDSpec(site_states=(("BD", 0), ("BA", 0)))
    esd_b = ESDSpec(site_states=(("BD", 0), ("BA", 0)))
    O = {("BD", 0, 0): 1.0}  # missing BA
    with pytest.raises(KeyError):
        esd_overlap(esd_a, esd_b, O)


def test_esd_overlap_different_fragments_raises():
    esd_a = ESDSpec(site_states=(("BD", 0),))
    esd_b = ESDSpec(site_states=(("BA", 0),))
    with pytest.raises(ValueError):
        esd_overlap(esd_a, esd_b, {})


# ------------------------------------------------------------ matrix

def test_esd_overlap_matrix_shape():
    esds_a = [
        ESDSpec(site_states=(("BD", 0), ("BA", 0))),
        ESDSpec(site_states=(("BD", 1), ("BA", 0))),
    ]
    esds_b = [
        ESDSpec(site_states=(("BD", 0), ("BA", 0))),
        ESDSpec(site_states=(("BD", 1), ("BA", 0))),
        ESDSpec(site_states=(("BD", 0), ("BA", 1))),
    ]
    O = {
        ("BD", 0, 0): 1.0, ("BD", 0, 1): 0.0,
        ("BD", 1, 0): 0.0, ("BD", 1, 1): 1.0,
        ("BA", 0, 0): 1.0, ("BA", 0, 1): 0.0,
    }
    M = esd_overlap_matrix(esds_a, esds_b, O)
    assert M.shape == (2, 3)


def test_esd_overlap_matrix_identity_case():
    """For an exact fragment identity, ESD overlap matrix is identity-like."""
    n = 3
    esds = [ESDSpec(site_states=(("BD", i), ("BA", 0))) for i in range(n)]
    # Fragment identity overlaps
    O = {}
    for i in range(n):
        for j in range(n):
            O[("BD", i, j)] = 1.0 if i == j else 0.0
    O[("BA", 0, 0)] = 1.0
    M = esd_overlap_matrix(esds, esds, O)
    np.testing.assert_allclose(M, np.eye(n), atol=1e-15)


# ------------------------------------------------------------ helpers

def test_fragment_overlaps_from_matrix():
    m = np.array([[0.9, 0.1], [0.05, 0.95]])
    O = fragment_overlaps_from_matrix("BD", m)
    assert O[("BD", 0, 0)] == 0.9
    assert O[("BD", 0, 1)] == 0.1
    assert O[("BD", 1, 0)] == 0.05
    assert O[("BD", 1, 1)] == 0.95


def test_merge_fragment_overlaps():
    d1 = {("BD", 0, 0): 1.0}
    d2 = {("BA", 0, 0): 0.9}
    merged = merge_fragment_overlaps(d1, d2)
    assert merged == {("BD", 0, 0): 1.0, ("BA", 0, 0): 0.9}


def test_complete_two_fragment_pipeline():
    """Simulate a two-fragment case end-to-end at the ESD level."""
    # Fragment overlap matrices (e.g. from wfoverlap at two geometries)
    O_BD = np.array([
        [0.999, 0.001, 0.000],
        [0.001, 0.998, 0.002],
        [0.000, 0.002, 0.997],
    ])
    O_BA = np.array([
        [0.998, 0.002],
        [0.002, 0.997],
    ])

    # Fragment overlaps dict
    frag_ov = merge_fragment_overlaps(
        fragment_overlaps_from_matrix("BD", O_BD),
        fragment_overlaps_from_matrix("BA", O_BA),
    )

    # ESDs: all combinations of (BD state, BA state)
    esds = [
        ESDSpec(site_states=(("BD", i), ("BA", j)))
        for i in range(3)
        for j in range(2)
    ]

    # The ESD overlap matrix should be the tensor product of the fragment matrices
    M = esd_overlap_matrix(esds, esds, frag_ov)
    assert M.shape == (6, 6)

    # Check tensor-product structure: M[(i,j),(i',j')] = O_BD[i,i'] * O_BA[j,j']
    for a, esd_a in enumerate(esds):
        for b, esd_b in enumerate(esds):
            i_a, j_a = esd_a.site_states[0][1], esd_a.site_states[1][1]
            i_b, j_b = esd_b.site_states[0][1], esd_b.site_states[1][1]
            expected = O_BD[i_a, i_b] * O_BA[j_a, j_b]
            assert abs(M[a, b] - expected) < 1e-15
