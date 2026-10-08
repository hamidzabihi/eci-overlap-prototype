"""
Phase correction for wavefunction overlap matrices.

The state-overlap matrix S_ij = <Psi_i^A | Psi_j^B> between two
geometries has arbitrary per-state phases: |Psi_i> and -|Psi_i>
describe the same physical state, but their overlap with a reference
state differs by a sign.

For NAMD, phases must be consistent across consecutive time steps so
that non-adiabatic couplings and hop probabilities are correctly
signed. This module implements two phase-correction algorithms:

1. **simple** - extract the phase of each diagonal element and divide
   it out. Fast and correct when states do not cross.

2. **robust** - minimize Tr(|log(U)|^2) with a Wigner-style sweep over
   pairs of states. Handles near-degenerate state pairs (trivial
   crossings) that break the simple method.

Both are ports of SHARC's ``utils.phase_correction`` (real) and
``utils.phase_correction_cmplx`` (complex) to a self-contained form
that we can test and reason about.

Reference
---------
Subotnik et al., J. Chem. Theory Comput. 2020, 16, 835-846.
"""

from __future__ import annotations
import numpy as np
from scipy import optimize


# ---------------------------------------------------------------- simple

def simple_phase_correction(S: np.ndarray) -> np.ndarray:
    """Return the per-state phase factors that make diag(S) real positive.

    Parameters
    ----------
    S : (n, n) complex overlap matrix

    Returns
    -------
    phases : (n,) complex array, |phases[i]| = 1
        Multiply column i of S by phases[i] to make S[i, i] real positive:
            S_corrected = S * phases[None, :]
        This convention matches the SHARC utility.
    """
    S = np.asarray(S)
    if S.ndim != 2 or S.shape[0] != S.shape[1]:
        raise ValueError("S must be square")
    n = S.shape[0]
    phases = np.ones(n, dtype=complex)
    for i in range(n):
        d = S[i, i]
        if abs(d) > 1e-12:
            phases[i] = np.conj(d) / abs(d)
    return phases


def apply_phases(S: np.ndarray, phases: np.ndarray) -> np.ndarray:
    """Multiply column i of S by phases[i]."""
    return S * phases[None, :]


# ---------------------------------------------------------------- robust

def _log_squared_norm(U: np.ndarray) -> float:
    """Compute Tr(|log(U)|^2) for a unitary matrix U."""
    eig = np.linalg.eigvals(U)
    return float(np.sum(np.abs(np.log(eig)) ** 2))


def robust_phase_correction(
    S: np.ndarray,
    *,
    tol: float = 1e-6,
    max_sweeps: int = 1000,
    seed: int | None = 0,
) -> np.ndarray:
    """Subotnik-style phase correction for complex overlap matrices.

    Iteratively minimizes Tr(|log(U)|^2) over pairs of phase factors
    applied directly to columns of S. The result is a set of phases
    that make the diagonal of S as real-positive as possible.

    Parameters
    ----------
    S : (n, n) complex matrix
    tol : convergence threshold on the pair-rotation angles
    max_sweeps : maximum number of full sweeps over pairs
    seed : RNG seed for the per-pair optimizer

    Returns
    -------
    phases : (n,) complex array, |phases[i]| = 1
    """
    S = np.asarray(S)
    if S.ndim != 2 or S.shape[0] != S.shape[1]:
        raise ValueError("S must be square")
    n = S.shape[0]

    S_work = S.copy()
    phases = np.ones(n, dtype=complex)

    # Initial per-column phase: make S_work[i, i] real and positive.
    for i in range(n):
        d = S_work[i, i]
        if abs(d) > 1e-12:
            ph = np.conj(d) / abs(d)
            S_work[:, i] *= ph
            phases[i] *= ph

    rng = np.random.default_rng(seed)

    def delta(theta, A1, A2, B1, B2):
        return (A1 * np.cos(theta) + A2 * np.cos(2 * theta)
                + B1 * np.sin(theta) + B2 * np.sin(2 * theta))

    sweeps = 0
    while True:
        done = True
        for j in range(n):
            for k in range(j + 1, n):
                U = S_work
                A1 = (6.0 * np.real(U[j, :] @ U[:, j] + U[k, :] @ U[:, k])
                      - 12.0 * np.real(U[j, k] * U[k, j])
                      - 6.0 * np.real(U[j, j] ** 2 + U[k, k] ** 2)
                      - 16.0 * np.real(U[j, j] + U[k, k]))
                A2 = 3.0 * np.real(U[j, j] ** 2 + U[k, k] ** 2)
                B1 = (6.0 * np.imag(U[k, :] @ U[:, k] - U[j, :] @ U[:, j])
                      - 6.0 * np.imag(U[k, k] ** 2 - U[j, j] ** 2)
                      - 16.0 * np.imag(U[k, k] - U[j, j]))
                B2 = 3.0 * np.imag(U[k, k] ** 2 - U[j, j] ** 2)

                seed_val = rng.uniform(0.0, 2.0 * np.pi)
                out = optimize.minimize(
                    delta, seed_val, args=(A1, A2, B1, B2),
                    method="BFGS", options={"gtol": 1e-8},
                )
                theta = float(out.x[0])

                S_work[:, j] *= np.exp(1j * theta)
                S_work[:, k] *= np.exp(-1j * theta)
                phases[j] *= np.exp(1j * theta)
                phases[k] *= np.exp(-1j * theta)

                if abs(theta) > tol:
                    done = False
        sweeps += 1
        if done or sweeps >= max_sweeps:
            break

    return phases


# ---------------------------------------------------------------- diagnostics

def max_diagonal_imag(S: np.ndarray) -> float:
    """Return the maximum |Im(S[i, i])|. A phase-corrected S has this ≈ 0."""
    return float(np.max(np.abs(S.diagonal().imag)))


def min_diagonal_magnitude(S: np.ndarray) -> float:
    """Return the minimum |S[i, i]|. A near-zero value warns of a trivial crossing."""
    return float(np.min(np.abs(S.diagonal())))
