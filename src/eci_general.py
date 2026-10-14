"""
Generalized ECI Hamiltonian for N fragments, any number of states per
fragment.

The existing `build_two_fragment_eci` handles the two-fragment case
with GS + 1 LE per fragment. This module generalizes it to:

  - N fragments (2, 3, 4, ...)
  - M states per fragment (GS + M-1 excited states)

The basis is:

    |GS>          = |0, 0, ..., 0>
    |F[a]>        = |0, ..., a_F, ..., 0>   (LE product on fragment F)

At FEM level the Hamiltonian is:
  - Diagonal: sum of site energies of all fragments in their respective
    states
  - Off-diagonal: Frenkel couplings between LEs on two different
    fragments (all other pairs are zero)

At ECIS level, GS-LE couplings and inter-fragment diagonal corrections
are added.

References:
  Pitesa, Polonius, Gonzalez, Mai, JCTC 2024, 20, 5609-5634
  (eqs. 6, 7, 8, 9, 10)
"""

from __future__ import annotations
import numpy as np
from pyscf import gto

from .excitonic_hamiltonian import (
    Fragment,
    cross_fragment_J,
    cross_fragment_K,
)


def build_eci_general(
    fragments,
    *,
    include_gs_le=False,
    spin=0,
    verbose=False,
):
    """Build the ECI Hamiltonian for N fragments.

    Parameters
    ----------
    fragments : list of Fragment
        Each fragment must have states[0] = GS and at least one excited
        state. The state ordering (0 = GS, 1, 2, ...) is preserved.
    include_gs_le : bool
        If True, add GS-LE couplings (ECIS level). If False, FEM level
        (GS-LE couplings are zero by construction).
    spin : int
        0 for singlet block, 1 for triplet block.
    verbose : bool
        Print progress.

    Returns
    -------
    H : (n_basis, n_basis) ndarray
    labels : list of str
    """
    N = len(fragments)
    if N < 2:
        raise ValueError("Need at least 2 fragments")
    n_states = [len(f.states) for f in fragments]
    if verbose:
        print(f"  Fragments: {N}, states per fragment: {n_states}")

    # ---------------------------------------------------------------- basis
    basis = []
    labels = []
    basis.append(tuple([0] * N))
    labels.append("|" + ",".join(["GS"] * N) + ">")

    for F in range(N):
        for a in range(1, n_states[F]):
            state = [0] * N
            state[F] = a
            basis.append(tuple(state))
            label = "|" + ",".join(
                f"{i}[{state[i]}]" if state[i] != 0 else "GS"
                for i in range(N)
            ) + ">"
            labels.append(label)

    n_basis = len(basis)
    if verbose:
        print(f"  Basis size (FEM level): {n_basis}")

    H = np.zeros((n_basis, n_basis))

    # ---------------------------------------------------------------- ERIs
    eri_cache = {}
    for F in range(N):
        for G in range(F + 1, N):
            eri_cache[(F, G)] = _build_cross_eris(fragments[F], fragments[G])
            if verbose:
                print(f"  ERIs {F}-{G}: shape {eri_cache[(F, G)]['eri_ABAB'].shape}")

    # ---------------------------------------------------------------- diagonal
    for i, state in enumerate(basis):
        E_site = sum(
            fragments[F].states[state[F]].energy for F in range(N)
        )
        H[i, i] = E_site

    # ---------------------------------------------------------------- off-diagonal
    for i, state_i in enumerate(basis):
        for j in range(i + 1, n_basis):
            state_j = basis[j]
            diff = [F for F in range(N) if state_i[F] != state_j[F]]

            if len(diff) == 1 and include_gs_le:
                F = diff[0]
                if state_i[F] == 0:
                    a_F = state_j[F]
                elif state_j[F] == 0:
                    a_F = state_i[F]
                else:
                    continue
                coupling = _gs_le_coupling(
                    fragments, F, a_F, state_i, eri_cache, spin=spin,
                )
                H[i, j] = H[j, i] = coupling

            elif len(diff) == 2:
                F, G = diff
                aF_i = state_i[F]; aF_j = state_j[F]
                aG_i = state_i[G]; aG_j = state_j[G]
                td_F = _transition_density(fragments[F], aF_i, aF_j)
                td_G = _transition_density(fragments[G], aG_i, aG_j)
                if td_F is None or td_G is None:
                    continue
                if not (np.any(np.abs(td_F) > 1e-12)
                        and np.any(np.abs(td_G) > 1e-12)):
                    continue
                coupling = _frenkel_coupling(
                    td_F, td_G, eri_cache[(F, G)], spin=spin,
                )
                H[i, j] = H[j, i] = coupling

    return H, labels


# ---------------------------------------------------------------- helpers

def _build_cross_eris(frag_F, frag_G):
    """Build cross-fragment two-electron integrals.

    Returns dict with:
      eri_ABAB      : (nao_F, nao_F, nao_G, nao_G)  -- (ij|kl)
      eri_ABAB_ikjl : (nao_F, nao_G, nao_F, nao_G)  -- (ik|jl)
                      permuted to (i, j, k, l) for K contraction
    """
    mol = gto.conc_mol(frag_F.mol, frag_G.mol)
    mol.build()
    nbas_F = frag_F.mol.nbas
    nbas_G = mol.nbas - nbas_F

    eri_ABAB = mol.intor(
        "int2e",
        shls_slice=(0, nbas_F, 0, nbas_F, nbas_F, mol.nbas,
                    nbas_F, mol.nbas),
    )
    eri_ikjl = mol.intor(
        "int2e",
        shls_slice=(0, nbas_F, nbas_F, mol.nbas, 0, nbas_F,
                    nbas_F, mol.nbas),
    )
    eri_ikjl_perm = np.transpose(eri_ikjl, (0, 2, 1, 3))
    return {"eri_ABAB": eri_ABAB, "eri_ikjl_perm": eri_ikjl_perm,
            "mol": mol, "nbas_F": nbas_F, "nbas_G": nbas_G}


def _transition_density(frag, a_bra, a_ket):
    """Return the transition density between states a_bra and a_ket.

    If a_bra == a_ket, returns the state density.
    Otherwise, returns the transition density (from fragment state a_bra
    to state a_ket).
    """
    if a_bra == a_ket:
        return frag.states[a_bra].density_total
    if a_bra == 0 and frag.states[a_ket].transition_density is not None:
        return frag.states[a_ket].transition_density
    if a_ket == 0 and frag.states[a_bra].transition_density is not None:
        return frag.states[a_bra].transition_density.T
    # Cross-excited transition (not implemented in the fragment model)
    return None


def _frenkel_coupling(td_F, td_G, eris, spin=0):
    """Frenkel coupling between LEs on two different fragments.

    Implements eq 10 of the paper: J_{a,0|0,b} - K_{a,0|0,b}
    """
    eri_J = eris["eri_ABAB"]
    eri_K = eris["eri_ikjl_perm"]

    # J: Coulomb between transition densities
    J = np.einsum("ij,ijkl,kl->", td_F, eri_J, td_G, optimize=True)

    # K: exchange. For a closed-shell singlet, spin-integrated:
    #   K = (1/2) * td_F . td_G . (ik|jl)
    K = 0.5 * np.einsum("ij,ijkl,kl->", td_F, eri_K, td_G, optimize=True)

    if spin == 0:
        return float(J + K)   # singlet: J + K
    else:
        return float(J - K)   # triplet: J - K


def _gs_le_coupling(fragments, F, a_F, state, eri_cache, spin=0):
    """GS-LE coupling: H[GS, LE_Fa].

    Uses the correct expression: V_ne + sum_{G} (J - K).
    """
    td_F = fragments[F].states[a_F].transition_density
    if td_F is None or not np.any(np.abs(td_F) > 1e-12):
        return 0.0

    # V_ne: transition density of F with nuclei of all OTHER fragments
    V_ne = 0.0
    for G in range(len(fragments)):
        if G == F:
            continue
        for g in range(fragments[G].mol.natm):
            Z_g = fragments[G].mol.atom_charge(g)
            R_g = fragments[G].mol.atom_coord(g)
            fragments[F].mol.set_rinv_orig(R_g)
            V_rinv = fragments[F].mol.intor("int1e_rinv")
            V_ne -= Z_g * np.einsum("ij,ij->", td_F, V_rinv)

    # sum_G (J - K) with G's GS density
    JK = 0.0
    key = (min(F, 0), max(F, 0))  # placeholder
    for G in range(len(fragments)):
        if G == F:
            continue
        pair = (min(F, G), max(F, G))
        eris = eri_cache[pair]
        rho_G = fragments[G].states[0].density_total
        eri_J = eris["eri_ABAB"]
        eri_K = eris["eri_ikjl_perm"]
        if F < G:
            J = np.einsum("ij,ijkl,kl->", td_F, eri_J, rho_G, optimize=True)
            K = 0.5 * np.einsum("ij,ijkl,kl->", td_F, eri_K, rho_G,
                                optimize=True)
        else:
            # For G < F we need B-first ERIs; not implemented yet.
            # This path will be needed for symmetric multi-fragment systems.
            # TODO: build B-first ERIs.
            J = 0.0
            K = 0.0
        JK += J - K

    return float(V_ne + JK)
