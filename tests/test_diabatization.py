"""
Tests for the state-assignment / diabatization module.
"""

import numpy as np
import pytest

from src.diabatization import (
    Assignment,
    assign_states,
    track_states,
    track_reference_state,
    assignment_table,
    character_matrix,
)


# ------------------------------------------------------------ assign_states

def test_assign_identity():
    """Diagonal overlap matrix -> trivial assignment."""
    S = np.eye(3, dtype=complex)
    a = assign_states(S)
    assert a.mapping == {0: 0, 1: 1, 2: 2}
    for j in range(3):
        assert a.overlaps[j] == 1.0


def test_assign_swap():
    """States 0 and 1 swapped between geometries."""
    S = np.array([
        [0.0, 1.0],
        [1.0, 0.0],
    ], dtype=complex)
    a = assign_states(S)
    assert a.mapping == {0: 1, 1: 0}


def test_assign_with_offdiagonal():
    """Clear maximum off the diagonal picks the right assignment."""
    S = np.array([
        [0.9, 0.1],
        [0.15, 0.95],
    ], dtype=complex)
    a = assign_states(S)
    assert a.mapping == {0: 0, 1: 1}


def test_assign_one_to_one():
    """Two new states both prefer reference 0 -> one is left unassigned."""
    S = np.array([
        [0.9, 0.8],
        [0.1, 0.05],
    ], dtype=complex)
    a = assign_states(S)
    # new state 0 wins (0.9 > 0.8), so new state 1 is unassigned
    assert a.mapping[0] == 0
    assert a.mapping[1] == -1
    assert a.unassigned_new == [1]


def test_assign_threshold():
    """Below-threshold overlaps are ignored."""
    S = np.array([
        [0.9, 0.001],
        [0.0, 0.0],
    ], dtype=complex)
    a = assign_states(S, threshold=0.1)
    assert a.mapping[0] == 0
    assert a.mapping[1] == -1


def test_assign_complex_phases():
    """Phases don't affect the assignment (uses |S|)."""
    S = np.array([
        [-0.9 + 0j, 0.05 + 0.02j],
        [0.05 - 0.02j, 0.85j],
    ], dtype=complex)
    a = assign_states(S)
    assert a.mapping == {0: 0, 1: 1}


def test_assign_more_new_than_ref():
    """3 new states, 2 references -> one new state unassigned."""
    S = np.array([
        [0.95, 0.1, 0.9],
        [0.05, 0.9, 0.1],
    ], dtype=complex)
    a = assign_states(S)
    assert len(a.reference_indices) <= 2
    assert len(a.unassigned_new) >= 1


# ------------------------------------------------------------ track_states

def test_track_states_length():
    S1 = np.eye(3, dtype=complex)
    S2 = np.eye(3, dtype=complex)
    a_list = track_states([S1, S2])
    assert len(a_list) == 2


def test_track_reference_state():
    """A reference state should show up at each step of the path."""
    S1 = np.eye(3, dtype=complex)
    S2 = np.array([
        [0.99, 0.01, 0.0],
        [0.01, 0.98, 0.02],
        [0.0, 0.02, 0.97],
    ], dtype=complex)
    a_list = track_states([S1, S2])
    track = track_reference_state(a_list, reference_index=0)
    assert track.reference_index == 0
    assert len(track.steps) == 2
    # step 0: matches new state 0 with overlap 1
    # step 1: matches new state 0 with overlap 0.99
    assert track.steps[0][1] == 0
    assert track.steps[1][1] == 0
    assert abs(track.steps[1][2] - 0.99) < 1e-12


# ------------------------------------------------------------ character_matrix

def test_character_matrix_shape():
    S1 = np.eye(3, dtype=complex)
    S2 = np.eye(3, dtype=complex)
    char = character_matrix([S1, S2], reference_index=0)
    assert char.shape == (2, 3)


def test_character_matrix_values():
    S1 = np.eye(2, dtype=complex)
    S2 = np.array([
        [0.9, 0.1],
        [0.1, 0.9],
    ], dtype=complex)
    char = character_matrix([S1, S2], reference_index=0)
    assert abs(char[0, 0] - 1.0) < 1e-12
    assert abs(char[0, 1] - 0.0) < 1e-12
    assert abs(char[1, 0] - 0.81) < 1e-12
    assert abs(char[1, 1] - 0.01) < 1e-12


# ------------------------------------------------------------ assignment_table

def test_assignment_table_readable():
    S = np.eye(3, dtype=complex)
    a = assign_states(S)
    table = assignment_table(a)
    assert "->" in table
    assert "S1" in table
    assert "S2" in table
    assert "S3" in table
