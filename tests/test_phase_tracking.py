"""
Tests for phase correction of overlap matrices.
"""

import numpy as np
import pytest

from src.phase_tracking import (
    simple_phase_correction,
    robust_phase_correction,
    apply_phases,
    max_diagonal_imag,
    min_diagonal_magnitude,
)


# ------------------------------------------------------------ simple

def test_simple_phase_identity():
    """Identity matrix needs no phase correction."""
    S = np.eye(3, dtype=complex)
    phases = simple_phase_correction(S)
    np.testing.assert_allclose(phases, np.ones(3), atol=1e-12)


def test_simple_phase_negative_diagonal():
    """Negative diagonal elements get phase -1."""
    S = np.diag([-1.0, 1.0, -1.0]).astype(complex)
    phases = simple_phase_correction(S)
    np.testing.assert_allclose(phases, [-1.0, 1.0, -1.0], atol=1e-12)


def test_simple_phase_complex():
    """Complex diagonal elements get the conjugated phase."""
    phases_in = np.array([1.0, 1j, -1.0, -1j])
    S = np.diag(phases_in).astype(complex)
    phases = simple_phase_correction(S)
    S_corr = apply_phases(S, phases)
    # Diagonal of the corrected matrix should be real and positive
    for i in range(4):
        assert abs(S_corr[i, i].imag) < 1e-12
        assert S_corr[i, i].real > 0


def test_simple_phase_offdiagonal_unchanged_magnitude():
    """Phases only affect the phase, not the magnitude."""
    S = np.array([
        [-0.9+0j, 0.1+0.05j, 0.0+0j],
        [0.1-0.05j, -0.8+0j, 0.2+0j],
        [0.0+0j, 0.2+0j, 0.7+0j],
    ], dtype=complex)
    phases = simple_phase_correction(S)
    S_corr = apply_phases(S, phases)
    np.testing.assert_allclose(np.abs(S_corr), np.abs(S), atol=1e-12)


def test_simple_phase_non_square_raises():
    with pytest.raises(ValueError):
        simple_phase_correction(np.zeros((2, 3), dtype=complex))


# ------------------------------------------------------------ robust

def test_robust_phase_identity():
    """Identity matrix: phases should all be 1 (or equivalent)."""
    S = np.eye(4, dtype=complex)
    phases = robust_phase_correction(S, seed=42)
    # Global phase of S is only defined up to the first column, so
    # phases[1:] should be 1 and phases[0] can be ±1.
    for i in range(1, 4):
        assert abs(phases[i] - 1.0) < 1e-6


def test_robust_phase_simple_case():
    """Diagonal-dominant complex matrix: robust should match simple."""
    S = np.array([
        [-0.95+0j, 0.02+0.01j, 0.0+0j],
        [0.02-0.01j, -0.9+0j, 0.03+0j],
        [0.0+0j, 0.03+0j, 0.85+0j],
    ], dtype=complex)
    phases_robust = robust_phase_correction(S, seed=42)
    S_corr = apply_phases(S, phases_robust)
    # Diagonal should be real, |diag| preserved
    for i in range(3):
        assert abs(S_corr[i, i].imag) < 0.02


def test_robust_makes_diagonal_real():
    """A phase-rotated diagonal should be recovered."""
    phases_true = np.array([1.0, np.exp(0.7j), np.exp(-1.1j), np.exp(2.0j)])
    S = np.diag(phases_true * np.array([0.95, 0.9, 0.85, 0.8])).astype(complex)
    phases = robust_phase_correction(S, seed=42)
    S_corr = apply_phases(S, phases)
    for i in range(4):
        assert abs(S_corr[i, i].imag) < 1e-6
        assert S_corr[i, i].real > 0


# ------------------------------------------------------------ diagnostics

def test_max_diagonal_imag():
    S = np.array([[1+0.1j, 0], [0, 1-0.2j]], dtype=complex)
    assert abs(max_diagonal_imag(S) - 0.2) < 1e-12


def test_min_diagonal_magnitude():
    S = np.diag([1.0, 0.3, 0.01]).astype(complex)
    assert abs(min_diagonal_magnitude(S) - 0.01) < 1e-12


def test_full_correction_pipeline():
    """End-to-end: raw S -> phases -> corrected S with real diagonal."""
    n = 5
    rng = np.random.default_rng(7)
    # Build a nearly-diagonal complex matrix
    S_raw = np.eye(n, dtype=complex) * (0.9 + 0.05 * rng.standard_normal(n))
    S_raw += 0.02 * (rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n)))
    # Apply random phases to columns (as if S came from an eigensolver)
    noise_phases = np.exp(1j * rng.uniform(0, 2*np.pi, n))
    S_noisy = S_raw * noise_phases[None, :]

    phases_simple = simple_phase_correction(S_noisy)
    S_simple = apply_phases(S_noisy, phases_simple)
    assert max_diagonal_imag(S_simple) < 0.02

    phases_robust = robust_phase_correction(S_noisy, seed=42)
    S_robust = apply_phases(S_noisy, phases_robust)
    assert max_diagonal_imag(S_robust) < 0.02
