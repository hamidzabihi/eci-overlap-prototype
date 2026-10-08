"""
Unit tests for the ECI state-overlap contraction.
"""

import numpy as np
import pytest

from src.state_overlap import (
    ecsf_overlap,
    state_overlap,
    state_overlap_simple,
    orthonormality_check,
    phase_correct,
)


# ------------------------------------------------------------ ecsf_overlap

def test_ecsf_overlap_identity_U():
    """With U = identity, O_ECSF = O_det."""
    n = 3
    U = np.eye(n)
    O_det = np.array([[0.9, 0.1, 0.0],
                      [0.1, 0.8, 0.2],
                      [0.0, 0.2, 0.7]])
    O_ECSF = ecsf_overlap(U, U, O_det)
    np.testing.assert_allclose(O_ECSF, O_det)


def test_ecsf_overlap_shape():
    U_A = np.eye(4)
    U_B = np.eye(4)
    O_det = np.eye(4)
    assert ecsf_overlap(U_A, U_B, O_det).shape == (4, 4)


def test_ecsf_overlap_incompatible_shapes_raises():
    with pytest.raises(ValueError):
        ecsf_overlap(np.eye(3), np.eye(3), np.eye(4))


# ------------------------------------------------------------ state_overlap

def test_state_overlap_identity_case():
    """If U = I, C = I, O_det = I, then S = I."""
    n = 3
    U = np.eye(n)
    C = np.eye(n)
    O_det = np.eye(n)
    S = state_overlap(U, C, U, C, O_det)
    np.testing.assert_allclose(S, np.eye(n), atol=1e-15)


def test_state_overlap_with_real_C():
    """A simple test with real coefficients."""
    U = np.eye(2)
    # State 0 = ECSF 0; state 1 = (ECSF 0 + ECSF 1)/sqrt(2)
    C = np.array([[1.0, 1/np.sqrt(2)],
                  [0.0, 1/np.sqrt(2)]])
    O_det = np.array([[1.0, 0.0],
                      [0.0, 1.0]])
    S = state_overlap(U, C, U, C, O_det)
    # S should equal C^T C
    np.testing.assert_allclose(S, C.T @ C, atol=1e-15)


def test_state_overlap_shape():
    nESD, nECSF_A, nECSF_B, nState_A, nState_B = 4, 3, 3, 2, 2
    U_A = np.random.default_rng(1).random((nESD, nECSF_A))
    U_B = np.random.default_rng(2).random((nESD, nECSF_B))
    C_A = np.random.default_rng(3).random((nECSF_A, nState_A))
    C_B = np.random.default_rng(4).random((nECSF_B, nState_B))
    O_det = np.random.default_rng(5).random((nESD, nESD))
    S = state_overlap(U_A, C_A, U_B, C_B, O_det)
    assert S.shape == (nState_A, nState_B)


def test_state_overlap_incompatible_shapes_raises():
    with pytest.raises(ValueError):
        state_overlap(np.eye(3), np.eye(2), np.eye(3), np.eye(2), np.eye(4))


# ------------------------------------------------------------ simple form

def test_state_overlap_simple_matches_full():
    """state_overlap_simple(UC_A, UC_B, O_det) == state_overlap(U_A, C_A, U_B, C_B, O_det)."""
    nESD, nECSF, nState = 4, 3, 2
    rng = np.random.default_rng(42)
    U_A = rng.random((nESD, nECSF))
    C_A = rng.random((nECSF, nState))
    U_B = rng.random((nESD, nECSF))
    C_B = rng.random((nECSF, nState))
    O_det = rng.random((nESD, nESD))

    S_full = state_overlap(U_A, C_A, U_B, C_B, O_det)
    S_simple = state_overlap_simple(U_A @ C_A, U_B @ C_B, O_det)
    np.testing.assert_allclose(S_full, S_simple, atol=1e-15)


# ------------------------------------------------------------ diagnostics

def test_orthonormality_check_identity():
    is_close, dev = orthonormality_check(np.eye(4))
    assert is_close
    assert dev < 1e-12


def test_orthonormality_check_non_identity():
    M = np.eye(3) + 0.1 * np.ones((3, 3))
    is_close, dev = orthonormality_check(M, tol=0.01)
    assert not is_close
    assert dev > 0.01


def test_orthonormality_check_non_square():
    is_close, dev = orthonormality_check(np.ones((3, 4)))
    assert not is_close


# ------------------------------------------------------------ phase correction

def test_phase_correct_makes_diagonal_real():
    S = np.array([
        [-0.9+0j, 0.1j],
        [0.1j,    -0.8+0j],
    ])
    S_corrected = phase_correct(S)
    # Diagonal should now be real and positive
    for i in range(2):
        assert abs(S_corrected[i, i].imag) < 1e-12
        assert S_corrected[i, i].real > 0
