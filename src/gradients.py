"""
Finite-difference gradients of the ECI Hamiltonian.

The ECI Hamiltonian depends on nuclear geometry through:

  - fragment site energies (RHF + CIS at each fragment geometry)
  - fragment MO coefficients (AO density and transition matrices
    used in the cross-fragment integrals)
  - cross-fragment two-electron integrals
  - nuclear attraction and repulsion terms

We compute the derivative of any eigenvalue with respect to any
Cartesian nuclear coordinate by central differences on the full
Hamiltonian. No analytic derivatives are used: at each perturbed
geometry the caller-provided function `rebuild_H` recomputes the
Hamiltonian from scratch.

This is a validation tool, not a production gradient. It costs
O(6 * N_atoms) Hamiltonian rebuilds per gradient, each of which
requires a fresh set of fragment SCF + CIS calculations. For
production work, analytic gradients would be needed; the derivative
terms for the EHF embedding charges and the fragment density
matrices are the non-trivial pieces.

The FD gradient is expected to agree with direct CIS gradients only
in direction, not in magnitude: the ECI and CIS wavefunctions are
different approximations to the same physics, and their potential
energy surfaces need not be parallel.
"""

from __future__ import annotations
from typing import Callable
import numpy as np


def fd_gradient(
    geom_a: list,
    geom_b: list,
    rebuild_H: Callable,
    *,
    delta: float = 5e-4,
    eigenvalue_index: int = 1,
) -> tuple[float, np.ndarray, np.ndarray]:
    """Central-difference gradient of one ECI eigenvalue.

    Parameters
    ----------
    geom_a, geom_b : list of [symbol, (x, y, z)]
        Reference geometries of the two fragments (Angstrom).
    rebuild_H : callable
        Function rebuild_H(geom_a, geom_b) -> H where H is the ECI
        Hamiltonian (n_basis, n_basis) in Hartree.
    delta : float
        Displacement for central differences (Angstrom).
    eigenvalue_index : int
        Index into the sorted eigenvalues. 0 = GS, 1 = S1, ...

    Returns
    -------
    e0 : float
        Reference eigenvalue (Hartree).
    grad_a : (N_a, 3) ndarray
        dE/dR in Hartree per Angstrom for each fragment-A atom.
    grad_b : (N_b, 3) ndarray
        Same for fragment B.
    """
    H0 = np.asarray(rebuild_H(geom_a, geom_b))
    e0 = float(np.linalg.eigvalsh(H0)[eigenvalue_index])

    grad_a = np.zeros((len(geom_a), 3))
    for i in range(len(geom_a)):
        for j in range(3):
            gp = _shift_atom(geom_a, i, j, +delta)
            gm = _shift_atom(geom_a, i, j, -delta)
            ep = float(np.linalg.eigvalsh(np.asarray(rebuild_H(gp, geom_b)))[eigenvalue_index])
            em = float(np.linalg.eigvalsh(np.asarray(rebuild_H(gm, geom_b)))[eigenvalue_index])
            grad_a[i, j] = (ep - em) / (2.0 * delta)

    grad_b = np.zeros((len(geom_b), 3))
    for i in range(len(geom_b)):
        for j in range(3):
            gp = _shift_atom(geom_b, i, j, +delta)
            gm = _shift_atom(geom_b, i, j, -delta)
            ep = float(np.linalg.eigvalsh(np.asarray(rebuild_H(geom_a, gp)))[eigenvalue_index])
            em = float(np.linalg.eigvalsh(np.asarray(rebuild_H(geom_a, gm)))[eigenvalue_index])
            grad_b[i, j] = (ep - em) / (2.0 * delta)

    return e0, grad_a, grad_b


def _shift_atom(geom, atom_idx, coord_idx, delta):
    """Return a shallow copy of geom with atom_idx shifted by delta."""
    out = []
    for k, (sym, xyz) in enumerate(geom):
        x, y, z = xyz
        if k == atom_idx:
            if coord_idx == 0:
                x = x + delta
            elif coord_idx == 1:
                y = y + delta
            elif coord_idx == 2:
                z = z + delta
        out.append([sym, (x, y, z)])
    return out
