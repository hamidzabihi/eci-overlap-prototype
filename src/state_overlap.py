"""
Compute ECI state overlaps from ESD overlaps.

Given:
  - O_det   : (nESD_A, nESD_B) matrix of ESD overlaps (from src.esd_overlap)
  - U_A, U_B: (nESD, nECSF) spin-adaptation matrices
  - C_A, C_B: (nECSF, nStates) ECI coefficient matrices

the ECSF overlap and state overlap matrices are:

    O_ECSF = U_A^T @ O_det @ U_B          (shape: nECSF_A x nECSF_B)
    S      = C_A^H @ O_ECSF @ C_B         (shape: nStates_A x nStates_B)

These are the quantities SHARC needs for surface-hopping NAMD (state
tracking, non-adiabatic couplings, hop probabilities).

The spin-adaptation matrix U and the ECI coefficients C come from the
ECI method itself (lib/ECI.py: ECIbasis[m].U and ECI.Psi['ECSF'][m]).
This module is agnostic to their origin: it only needs the matrices.
"""

from __future__ import annotations
import numpy as np


def ecsf_overlap(
    U_A: np.ndarray,
    U_B: np.ndarray,
    O_det: np.ndarray,
) -> np.ndarray:
    """Contract ESD overlaps with spin-adaptation matrices.

    Parameters
    ----------
    U_A   : (nESD_A, nECSF_A)  spin-adaptation matrix at geometry A
    U_B   : (nESD_B, nECSF_B)  spin-adaptation matrix at geometry B
    O_det : (nESD_A, nESD_B)   ESD overlap matrix

    Returns
    -------
    O_ECSF : (nECSF_A, nECSF_B)  ECSF overlap matrix
    """
    U_A = np.asarray(U_A)
    U_B = np.asarray(U_B)
    O_det = np.asarray(O_det)
    if U_A.shape[0] != O_det.shape[0]:
        raise ValueError(
            f"U_A has {U_A.shape[0]} ESDs, O_det has {O_det.shape[0]} rows"
        )
    if U_B.shape[0] != O_det.shape[1]:
        raise ValueError(
            f"U_B has {U_B.shape[0]} ESDs, O_det has {O_det.shape[1]} cols"
        )
    return U_A.T @ O_det @ U_B


def state_overlap(
    U_A: np.ndarray,
    C_A: np.ndarray,
    U_B: np.ndarray,
    C_B: np.ndarray,
    O_det: np.ndarray,
) -> np.ndarray:
    """Compute the full state-overlap matrix S_ij = <Psi_i^A | Psi_j^B>.

    Parameters
    ----------
    U_A, U_B : (nESD, nECSF)   spin-adaptation matrices
    C_A, C_B : (nECSF, nState) ECI coefficient matrices
    O_det    : (nESD, nESD)    ESD overlap matrix

    Returns
    -------
    S : (nState_A, nState_B) complex matrix
    """
    O_ECSF = ecsf_overlap(U_A, U_B, O_det)
    C_A = np.asarray(C_A)
    C_B = np.asarray(C_B)
    if C_A.shape[0] != O_ECSF.shape[0]:
        raise ValueError(
            f"C_A has {C_A.shape[0]} ECSFs, O_ECSF has {O_ECSF.shape[0]} rows"
        )
    if C_B.shape[0] != O_ECSF.shape[1]:
        raise ValueError(
            f"C_B has {C_B.shape[0]} ECSFs, O_ECSF has {O_ECSF.shape[1]} cols"
        )
    return C_A.conj().T @ O_ECSF @ C_B


# ---------------------------------------------------------------- helpers

def state_overlap_simple(
    UC_A: np.ndarray,
    UC_B: np.ndarray,
    O_det: np.ndarray,
) -> np.ndarray:
    """Same as state_overlap but taking the pre-multiplied (U @ C) matrices.

    This is useful when U and C are not separately available, or when the
    ECI code stores the combined transformation directly. The result is
    identical to state_overlap(U_A, C_A, U_B, C_B, O_det) when
    UC_A = U_A @ C_A and UC_B = U_B @ C_B.

    Parameters
    ----------
    UC_A, UC_B : (nESD, nState)  combined transformation matrices
    O_det      : (nESD, nESD)   ESD overlap matrix

    Returns
    -------
    S : (nState_A, nState_B)
    """
    UC_A = np.asarray(UC_A)
    UC_B = np.asarray(UC_B)
    O_det = np.asarray(O_det)
    return UC_A.conj().T @ O_det @ UC_B


def orthonormality_check(S: np.ndarray, tol: float = 1e-8) -> tuple[bool, float]:
    """Check whether a state-overlap matrix is close to the identity.

    Returns
    -------
    (is_close, max_deviation)
    """
    if S.shape[0] != S.shape[1]:
        return False, float("inf")
    dev = np.max(np.abs(S - np.eye(S.shape[0])))
    return dev < tol, dev


def phase_correct(S: np.ndarray) -> np.ndarray:
    """Apply a simple phase correction so that the diagonal of S is real positive.

    The overlap matrix between two adiabatic states at nearby geometries
    should have unit-modulus diagonal entries when the states match.
    Multiply each column by the conjugate phase of its diagonal element
    to make the diagonal real. This is the "simple" phase convention
    (see utils.phase_correction_cmplx in SHARC).

    Note: only meaningful when S[i,i] is non-negligible.
    """
    S = np.array(S, dtype=complex)
    n = S.shape[0]
    for i in range(n):
        d = S[i, i]
        if abs(d) > 1e-12:
            S[:, i] *= np.conj(d) / abs(d)
    return S
