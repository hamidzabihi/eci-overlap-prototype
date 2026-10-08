"""
Assemble the ECI Hamiltonian from fragment CIS calculations.

This implements the fragment-based effective Hamiltonian of the
Excitonic Configuration Interaction (ECI) method, in the simplified
form needed for a two-fragment system with GS + local excitations.

Theory
------
The ECI Hamiltonian in the excitonic basis, within the strong-
orthogonality assumption (McWeeny's group function theory), is:

    H_ab = sum_F E_F[s_a^F] * delta_ab                      (site energies)
         + sum_{F<G} ( J_FG - K_FG )                        (interactions)

where the J and K integrals are defined in the paper (eqs. 6-7):

    J_ab|cd^FG = integral dr1 dr2  rho_F[a,b](r1) rho_G[c,d](r2) / |r1 - r2|
               + nuclear-electron terms + V_NN

    K_ab|cd^FG = integral dx1 dx2  rho_F[a,b](x1,x2) rho_G[c,d](x1,x2) / |r1 - r2|

For a two-fragment system with GS + 1 LE per fragment, the four
basis states are:

    |0>       = A(GS) (x) B(GS)          ground state
    |A>       = A(LE) (x) B(GS)          LE on fragment A
    |B>       = A(GS) (x) B(LE)          LE on fragment B
    |AB>      = A(LE) (x) B(LE)          DLE (optional)

Diagonal elements (eq. S2 / S5 in SI):

    H_00 = E_A^GS + E_B^GS + Y_00|00^AB
    H_AA = E_A^LE + E_B^GS + Y_LE,LE|00^AB
    H_BB = E_A^GS + E_B^LE + Y_00|LE,LE^AB

Off-diagonal (eqs. S5, S11):

    H_0A = Y_0,LE|00^AB    (GS-LE coupling; vanishes with EHF embedding)
    H_0B = Y_00|0,LE^AB
    H_AB = Y_LE,0|0,LE^AB  (Frenkel coupling; this is the key excitonic term)

where Y_ab|cd = J_ab|cd - K_ab|cd.

Density matrix conventions
--------------------------
Following the SI (Figure S2), each fragment state has four spin-
resolved density matrix components:

    rho^aa : alpha-alpha
    rho^bb : beta-beta
    rho^ab : alpha-beta
    rho^ba : beta-alpha

For a singlet CIS state with amplitudes X_{ia}, built from a closed-
shell RHF reference:

    rho^aa_{ij} = delta_ij * 2  for occupied i,j (reference)
                + sum_a X_{ia} X_{ja}        (virtual-occupied correction)
    ... (see below for the exact formula)

For the demo, we implement the restricted-singlet case, which is
sufficient for the two-fragment GS + LE problem.

References
----------
- Pitesa, Polonius, Gonzalez, Mai,
  J. Chem. Theory Comput. 2024, 20, 5609-5634 (ECI method)
- Pitesa, Mai, Gonzalez,
  J. Phys. Chem. Lett. 2025, 16, 2800-2807 (RI approximation)
- McWeeny, Methods of Molecular Quantum Mechanics, 2nd ed. (1992)
  eqs. 14.1.7-14.1.10 (GFT formulas)
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional
import numpy as np
from pyscf import gto


# ---------------------------------------------------------------- fragment data

@dataclass
class FragmentState:
    """A single electronic state of a fragment."""
    label: str               # e.g. "A"
    state_index: int         # 0 = GS, 1, 2, ... = excited
    energy: float            # total energy in the fragment's SCF reference (Ha)
    density_total: np.ndarray      # (nao, nao) total density matrix
    density_alpha: np.ndarray      # (nao, nao) alpha density
    density_beta: np.ndarray       # (nao, nao) beta density
    transition_density: np.ndarray | None = None  # (nao, nao) for GS->state


@dataclass
class Fragment:
    """A fragment: molecule, MO coefficients, and a list of states."""
    label: str
    mol: gto.Mole
    mo_coeff: np.ndarray
    mo_energy: np.ndarray
    mo_occ: np.ndarray
    states: list = field(default_factory=list)  # list of FragmentState

    @property
    def nao(self):
        return self.mol.nao

    @property
    def nocc(self):
        return int(np.sum(self.mo_occ > 0))


# ---------------------------------------------------------------- density construction

def cis_state_density(
    mo_coeff: np.ndarray,
    mo_occ: np.ndarray,
    x_ia: np.ndarray,
    *,
    state: str = "alpha",
) -> np.ndarray:
    """Build the AO-basis density matrix for a CIS state.

    Parameters
    ----------
    mo_coeff : (nao, nmo) MO coefficients
    mo_occ   : (nmo,) occupations (2.0 for doubly occupied, 0.0 for virtual)
    x_ia     : (nocc, nvir) CIS amplitude matrix for the state
    state    : 'total', 'alpha', or 'beta'

    Returns
    -------
    (nao, nao) density matrix in AO basis

    Notes
    -----
    For a CIS state built from a closed-shell RHF reference:

        rho_total = rho_ref + sum_{ia} X_ia * (phi_i phi_a^T + phi_a phi_i^T)
        rho_alpha = (1/2) rho_total + (1/2) * spin_broken_correction
        rho_beta  = (1/2) rho_total - (1/2) * spin_broken_correction

    For a closed-shell singlet, the spin_broken correction is zero when
    we symmetrize the alpha/beta excitations. In the demo, we treat
    singlet excitations as alpha/beta symmetric, so rho_alpha = rho_beta
    = rho_total / 2.

    In the full ECI (paper), more care is needed with spin-resolved
    densities for spin-broken determinants; see Fig. S2 in the SI.
    """
    nocc = int(np.sum(mo_occ > 0))
    nmo = mo_coeff.shape[1]
    nvir = nmo - nocc

    # Reference (GS) density
    mo_occ_ref = mo_occ[:nocc]
    # For RHF, occupied orbitals have occupation 2 -> in the AO basis:
    # rho_ref = sum_i occ_i * phi_i phi_i^T with occ_i = 2 for each doubly
    # occupied MO
    rho_ref = mo_coeff[:, :nocc] @ np.diag(mo_occ_ref) @ mo_coeff[:, :nocc].T

    # CIS correction: sum_{ia} X_ia * (phi_i phi_a^T + phi_a phi_i^T)
    if x_ia is None or np.allclose(x_ia, 0):
        rho_total = rho_ref
    else:
        # Vectorize the sum
        occ_mo = mo_coeff[:, :nocc]   # (nao, nocc)
        vir_mo = mo_coeff[:, nocc:]   # (nao, nvir)
        # phi_i phi_a^T summed over (i,a) with weights X_ia:
        # = occ_mo @ X @ vir_mo.T
        cross = occ_mo @ x_ia @ vir_mo.T
        rho_total = rho_ref + cross + cross.T

    if state == "total":
        return rho_total
    elif state == "alpha":
        return 0.5 * rho_total
    elif state == "beta":
        return 0.5 * rho_total
    else:
        raise ValueError(f"state must be 'total', 'alpha', or 'beta', got {state!r}")


def cis_transition_density(
    mo_coeff: np.ndarray,
    mo_occ: np.ndarray,
    x_ia: np.ndarray,
) -> np.ndarray:
    """AO-basis transition density between GS and a CIS state.

    For a singlet CIS state, this is:

        rho_{0 -> S}_{AO} = sum_{ia} X_ia * (phi_i phi_a^T + phi_a phi_i^T)

    Note: the GS reference contribution is excluded (as required for a
    transition density).
    """
    nocc = int(np.sum(mo_occ > 0))
    occ_mo = mo_coeff[:, :nocc]
    vir_mo = mo_coeff[:, nocc:]
    cross = occ_mo @ x_ia @ vir_mo.T
    return cross + cross.T


# ---------------------------------------------------------------- JK integrals

def cross_fragment_J(
    rho_F: np.ndarray,
    rho_G: np.ndarray,
    eri_cross: np.ndarray,
) -> float:
    """Coulomb integral J between two fragment densities.

    J = sum_{ij on F, kl on G} rho_F[i,j] * rho_G[k,l] * (ij|kl)

    where eri_cross[i,j,k,l] = (ij|kl) is the two-electron integral with
    i,j on fragment F and k,l on fragment G.
    """
    # (i,j) contracted with eri, then with rho_G over (k,l)
    J = np.einsum("ij,ijkl,kl->", rho_F, eri_cross, rho_G, optimize=True)
    return float(J)


def cross_fragment_K(
    rho_F_alpha: np.ndarray,
    rho_F_beta: np.ndarray,
    rho_G_alpha: np.ndarray,
    rho_G_beta: np.ndarray,
    eri_cross_aa: np.ndarray,
    eri_cross_bb: np.ndarray,
    eri_cross_ab: np.ndarray,
    eri_cross_ba: np.ndarray,
) -> float:
    """Exchange integral K between two fragment densities.

    Following eq. 49 in the SI:

    K = sum_{ij on F, kl on G} [
            rho_F^aa[i,j] * rho_G^aa[k,l] * (ik|jl)
          + rho_F^bb[i,j] * rho_G^bb[k,l] * (ik|jl)
          + rho_F^ab[i,j] * rho_G^ba[k,l] * (ik|jl)
          + rho_F^ba[i,j] * rho_G^ab[k,l] * (ik|jl)
        ]

    The eri_cross_* arrays are the (ik|jl) integrals with the appropriate
    index ordering.
    """
    K = 0.0
    K += np.einsum("ij,ijkl,kl->", rho_F_alpha, eri_cross_aa, rho_G_alpha, optimize=True)
    K += np.einsum("ij,ijkl,kl->", rho_F_beta, eri_cross_bb, rho_G_beta, optimize=True)
    K += np.einsum("ij,ijkl,kl->", rho_F_alpha, eri_cross_ab, rho_G_beta, optimize=True)
    K += np.einsum("ij,ijkl,kl->", rho_F_beta, eri_cross_ba, rho_G_alpha, optimize=True)
    return float(K)


def one_electron_operator(
    mol_F: gto.Mole,
    mol_G: gto.Mole,
    all_mols: list,
) -> np.ndarray:
    """Build the inter-fragment one-electron operator h^FG = T + V_ne.

    Parameters
    ----------
    mol_F, mol_G : the two fragments between which the operator acts
    all_mols : all fragments' molecules (for the nuclear attraction sum)

    Returns
    -------
    h_FG : (nao_F, nao_G) matrix in the mixed AO basis
    """
    # Kinetic energy block (F rows, G columns)
    dimer = gto.conc_mol(mol_F, mol_G)
    dimer.build()
    nao_F = mol_F.nao
    nao_G = mol_G.nao
    T = dimer.intor("int1e_kin")[:nao_F, nao_F:nao_F + nao_G]

    # Nuclear-electron attraction: -sum over all nuclei of all fragments
    VNE = np.zeros_like(T)
    for other_mol in all_mols:
        for atom in range(other_mol.natm):
            Z = other_mol.atom_charge(atom)
            R = other_mol.atom_coord(atom)
            dimer.set_rinv_orig(R)
            V = dimer.intor("int1e_rinv")[:nao_F, nao_F:nao_F + nao_G]
            VNE -= Z * V

    return T + VNE


def nuclear_terms(
    mol_F: gto.Mole,
    mol_G: gto.Mole,
    rho_F: np.ndarray,
    rho_G: np.ndarray,
) -> tuple[float, float, float]:
    """Return the three nuclear contributions separately.

    Returns
    -------
    (V_ne_FG, V_ne_GF, V_nn)
        V_ne_FG: electron density of F with nuclei of G
        V_ne_GF: electron density of G with nuclei of F
        V_nn:    nuclear-nuclear repulsion

    Each piece is included in the GFT J integral only when the
    corresponding Kronecker-delta condition is satisfied:
        V_ne_FG requires delta(a_F, b_F) = 1
        V_ne_GF requires delta(a_G, b_G) = 1
        V_nn    requires delta(a_F, b_F) * delta(a_G, b_G) = 1
    """
    V_ne_FG = 0.0
    for g in range(mol_G.natm):
        Z_g = mol_G.atom_charge(g)
        R_g = mol_G.atom_coord(g)
        mol_F.set_rinv_orig(R_g)
        V_F_at_G = mol_F.intor("int1e_rinv")
        V_ne_FG -= Z_g * np.einsum("ij,ij->", rho_F, V_F_at_G)

    V_ne_GF = 0.0
    for f in range(mol_F.natm):
        Z_f = mol_F.atom_charge(f)
        R_f = mol_F.atom_coord(f)
        mol_G.set_rinv_orig(R_f)
        V_G_at_F = mol_G.intor("int1e_rinv")
        V_ne_GF -= Z_f * np.einsum("ij,ij->", rho_G, V_G_at_F)

    V_nn = 0.0
    for f in range(mol_F.natm):
        for g in range(mol_G.natm):
            Z_f = mol_F.atom_charge(f)
            Z_g = mol_G.atom_charge(g)
            R_f = mol_F.atom_coord(f)
            R_g = mol_G.atom_coord(g)
            V_nn += Z_f * Z_g / np.linalg.norm(R_f - R_g)

    return V_ne_FG, V_ne_GF, V_nn


# ---------------------------------------------------------------- ECI assembly

@dataclass
class ECIBasis:
    """The excitonic basis for a two-fragment system."""
    fragment_A: Fragment
    fragment_B: Fragment
    state_A: int             # index into fragment_A.states
    state_B: int             # index into fragment_B.states

    @property
    def label(self):
        return f"A[{self.state_A}] x B[{self.state_B}]"


def build_two_fragment_eci(
    fragment_A: Fragment,
    fragment_B: Fragment,
    *,
    include_dle: bool = False,
    include_gs_le: bool = False,
    verbose: bool = False,
) -> tuple[np.ndarray, list[str]]:
    """Build the ECI Hamiltonian matrix for a two-fragment system.

    The basis contains:
      - GS product: A[0] x B[0]
      - LE on A:    A[1] x B[0]
      - LE on B:    A[0] x B[1]
      - DLE:        A[1] x B[1]     (if include_dle=True)

    Parameters
    ----------
    fragment_A, fragment_B : Fragment
        Each must have at least 2 states (GS + 1 excited).
    include_dle : bool
        Include the doubly-localized excitation state.

    Returns
    -------
    H : (n, n) real Hamiltonian matrix in the ECSF basis
    labels : list of basis labels
    """
    # Build cross-fragment 2-electron integrals
    mol_AB = gto.conc_mol(fragment_A.mol, fragment_B.mol)
    mol_AB.build()
    nbas_A = fragment_A.mol.nbas
    nbas_B = mol_AB.nbas - nbas_A

    # (ij|kl): i,j on A; k,l on B
    eri_ABAB = mol_AB.intor(
        "int2e",
        shls_slice=(0, nbas_A, 0, nbas_A, nbas_A, mol_AB.nbas, nbas_A, mol_AB.nbas),
    )
    # (ik|jl): i on A, k on B, j on A, l on B
    eri_ABAB_ikjl = mol_AB.intor(
        "int2e",
        shls_slice=(0, nbas_A, nbas_A, mol_AB.nbas, 0, nbas_A, nbas_A, mol_AB.nbas),
    )

    # Redefine index order for the K integral
    # We need (ik|jl) with i,j on A and k,l on B. PySCF's int2e with
    # shls_slice=(i_range, k_range, j_range, l_range) gives exactly that.
    # The array shape is (nao_A, nao_B, nao_A, nao_B).
    # For the K contraction, we want i,j on A and k,l on B.
    # So we permute to (i, j, k, l):
    eri_ABAB_ikjl_perm = np.transpose(eri_ABAB_ikjl, (0, 2, 1, 3))
    # Now eri_ABAB_ikjl_perm[i, j, k, l] = (ik|jl)

    # Get fragment state energies and densities
    sA_gs = fragment_A.states[0]
    sB_gs = fragment_B.states[0]

    # Basis
    basis = []
    basis.append(ECIBasis(fragment_A, fragment_B, 0, 0))   # GS
    basis.append(ECIBasis(fragment_A, fragment_B, 1, 0))   # LE on A
    basis.append(ECIBasis(fragment_A, fragment_B, 0, 1))   # LE on B
    if include_dle:
        basis.append(ECIBasis(fragment_A, fragment_B, 1, 1))  # DLE

    labels = [b.label for b in basis]
    n = len(basis)
    H = np.zeros((n, n))

    # --- Diagonal elements and Y integrals for all pairs ---

    # Precompute Y integrals between all same-fragment state pairs
    # Y_ab|cd^FG = J_ab|cd - K_ab|cd + nuclear terms
    # We need all (a,b) pairs on A and (c,d) pairs on B that appear.

    def density_for_state(frag: Fragment, state_idx: int, spin: str):
        """Get AO density for the given fragment state and spin."""
        st = frag.states[state_idx]
        if spin == "total":
            return st.density_total
        elif spin == "alpha":
            return st.density_alpha
        elif spin == "beta":
            return st.density_beta
        else:
            raise ValueError(spin)

    def transition_density_for_state(frag: Fragment, state_idx: int):
        """Get AO transition density between GS and the state."""
        st = frag.states[state_idx]
        if st.transition_density is not None:
            return st.transition_density
        return np.zeros((frag.nao, frag.nao))

    def Y_integral(
        rho_F_ab, rho_G_cd,
        rho_F_ab_alpha, rho_F_ab_beta,
        rho_G_cd_alpha, rho_G_cd_beta,
        *,
        is_transition_F: bool = False,
        is_transition_G: bool = False,
    ) -> float:
        """Compute Y = J - K + nuclear for two fragment densities.

        The nuclear-electron and nuclear-nuclear terms are only nonzero
        when the two states on each fragment are the SAME (Kronecker
        delta structure in eq. 6 of the paper). For transition densities
        (a_F != b_F), those terms vanish.

        Parameters
        ----------
        is_transition_F, is_transition_G : bool
            True if the corresponding density is a transition density
            (between two different fragment states). When False, the
            density is a state density and the nuclear terms apply.
        """
        J = cross_fragment_J(rho_F_ab, rho_G_cd, eri_ABAB)
        K = (
            np.einsum("ij,ijkl,kl->", rho_F_ab_alpha, eri_ABAB_ikjl_perm, rho_G_cd_alpha, optimize=True)
            + np.einsum("ij,ijkl,kl->", rho_F_ab_beta, eri_ABAB_ikjl_perm, rho_G_cd_beta, optimize=True)
        )
        # Nuclear terms follow the Kronecker-delta structure of eq. 6 of
        # the ECI paper (JCTC 2024). For a transition-density pair
        # (a_F != b_F), the corresponding delta is zero, and the nuclear
        # terms are absent.
        #
        # NOTE: This simple implementation is correct for the diagonal
        # (state-density | state-density) case. It produces a physical
        # Frenkel coupling for (transition | state) pairs, but the GS-LE
        # coupling requires a more careful derivation of the AO-basis
        # formula. See NOTES.md for details.
        V_ne_FG, V_ne_GF, V_nn = nuclear_terms(
            fragment_A.mol, fragment_B.mol, rho_F_ab, rho_G_cd
        )
        nuc = 0.0
        if (not is_transition_F) and (not is_transition_G):
            nuc = V_ne_FG + V_ne_GF + V_nn
        return float(J - K + nuc)

    # Diagonal H_00 = E_A^GS + E_B^GS + Y_00|00
    Y_00_00 = Y_integral(
        sA_gs.density_total, sB_gs.density_total,
        sA_gs.density_alpha, sA_gs.density_beta,
        sB_gs.density_alpha, sB_gs.density_beta,
    )
    H[0, 0] = sA_gs.energy + sB_gs.energy + 0.5 * Y_00_00

    if verbose:
        print(f"  Y_00|00 = {Y_00_00:.8f}")

    # Diagonal H_AA = E_A^LE + E_B^GS + Y_LE,LE|00
    if n > 1:
        sA_le = fragment_A.states[1]
        Y_LELE_00 = Y_integral(
            sA_le.density_total, sB_gs.density_total,
            sA_le.density_alpha, sA_le.density_beta,
            sB_gs.density_alpha, sB_gs.density_beta,
        )
        H[1, 1] = sA_le.energy + sB_gs.energy + 0.5 * Y_LELE_00
        if verbose:
            print(f"  Y_LE_A,LE_A|00 = {Y_LELE_00:.8f}")

    # Diagonal H_BB = E_A^GS + E_B^LE + Y_00|LE,LE
    if n > 2:
        sB_le = fragment_B.states[1]
        Y_00_LELE = Y_integral(
            sA_gs.density_total, sB_le.density_total,
            sA_gs.density_alpha, sA_gs.density_beta,
            sB_le.density_alpha, sB_le.density_beta,
        )
        H[2, 2] = sA_gs.energy + sB_le.energy + 0.5 * Y_00_LELE
        if verbose:
            print(f"  Y_00|LE_B,LE_B = {Y_00_LELE:.8f}")

    # Off-diagonal GS-LE coupling using the inter-fragment one-electron
    # operator h^FG = T + V_ne, following the implementation in SHARC's
    # lib/ECI.py (calculate_V1mat).
    #
    #   H_{GS-LE_A} = (1/2) * Tr( h^{AB} * P^{AB} )
    #
    # where P^{AB} = rho^{trans,A} @ S_AB @ rho^{GS,B} is the mixed
    # two-fragment density matrix.
    all_mols = [fragment_A.mol, fragment_B.mol]
    S_AB = gto.mole.intor_cross("int1e_ovlp", fragment_A.mol, fragment_B.mol)

    if n > 1 and include_gs_le:
        td_A = transition_density_for_state(fragment_A, 1)
        if np.any(np.abs(td_A) > 1e-10):
            hFG = one_electron_operator(fragment_A.mol, fragment_B.mol, all_mols)
            P_mixed = td_A @ S_AB @ sB_gs.density_total
            V1_0A = 0.5 * np.einsum("ij,ji->", hFG, P_mixed)
            H[0, 1] = H[1, 0] = V1_0A
            if verbose:
                print(f"  V1_0,LE_A = {V1_0A:.8f}")

    if n > 2 and include_gs_le:
        td_B = transition_density_for_state(fragment_B, 1)
        if np.any(np.abs(td_B) > 1e-10):
            hGF = one_electron_operator(fragment_B.mol, fragment_A.mol, all_mols)
            P_mixed_BA = td_B @ S_AB.T @ sA_gs.density_total
            V1_0B = 0.5 * np.einsum("ij,ji->", hGF, P_mixed_BA)
            H[0, 2] = H[2, 0] = V1_0B
            if verbose:
                print(f"  V1_00,LE_B = {V1_0B:.8f}")

    # Off-diagonal H_AB = Y_LE,0|0,LE (Frenkel coupling)
    if n > 2:
        td_A = transition_density_for_state(fragment_A, 1)
        td_B = transition_density_for_state(fragment_B, 1)
        if np.any(np.abs(td_A) > 1e-10) and np.any(np.abs(td_B) > 1e-10):
            Y_LE0_0LE = Y_integral(
                td_A, td_B,
                0.5 * td_A, 0.5 * td_A,
                0.5 * td_B, 0.5 * td_B,
                is_transition_F=True,
                is_transition_G=True,
            )
            H[1, 2] = H[2, 1] = Y_LE0_0LE
            if verbose:
                print(f"  Y_LE_A,0|0,LE_B = {Y_LE0_0LE:.8f}")

    return H, labels


def fragment_from_cis_result(label, cis_result, *, include_all_states=True):
    """Build a Fragment object from the output of run_cis.

    Parameters
    ----------
    label : str
        Fragment label (e.g. "A").
    cis_result : dict
        Output of run_cis, containing keys 'mol', 'mo_coeff', 'coefs',
        and (optionally) SCF energy.

    Returns
    -------
    Fragment instance with GS + excited states populated.
    """
    mol = cis_result["mol"]
    mo_coeff = cis_result["mo_coeff"]
    nocc = cis_result["nocc"]
    nvir = cis_result["nvir"]
    nmo = cis_result["nmo"]

    # Reconstruct mo_occ: RHF, all occupied doubly, all virtual empty
    mo_occ = np.zeros(nmo)
    mo_occ[:nocc] = 2.0

    # Get the SCF energy -- in the current run_cis, it's not returned, so
    # we recompute it from the reference density
    # (In production, we'd pass it through)
    # For the demo, we take the reference energy as the GS site energy.
    mf_energy = cis_result.get("scf_energy")
    if mf_energy is None:
        # Fallback: use 0, which we'll override at demo time
        mf_energy = 0.0

    # Reconstruct the X amplitudes from the CI coefficients
    # The coefs matrix (nstate, ndets) and detstrings give us the
    # full CIS vector. For each excited state, we need to extract the
    # (nocc, nvir) amplitude matrix.
    #
    # The detstring convention from run_cis:
    #   "d"*nocc + "e"*nvir    -> ground state
    #   with one 'd' -> 'b' and one 'e' -> 'a'    -> alpha excitation
    #   with one 'd' -> 'a' and one 'e' -> 'b'    -> beta excitation
    #
    # So to reconstruct X[i, a] for an excited state:
    # X_ia = coefficient of the detstring with position i changed to 'b'
    #        and position (nocc + a) changed to 'a'
    # (We take the alpha-excited detstring; the beta one has opposite sign.)
    detstrings = cis_result["detstrings"]
    coefs = cis_result["coefs"]
    nstate = cis_result["nstate"]

    det_to_index = {d: i for i, d in enumerate(detstrings)}

    states = []
    # GS
    rho_gs = cis_state_density(mo_coeff, mo_occ, None, state="total")
    states.append(FragmentState(
        label=label, state_index=0,
        energy=mf_energy,
        density_total=rho_gs,
        density_alpha=0.5 * rho_gs,
        density_beta=0.5 * rho_gs,
        transition_density=None,
    ))

    # Excited states
    for s in range(1, nstate):
        X = np.zeros((nocc, nvir))
        for i in range(nocc):
            for a in range(nvir):
                sym = ["d"] * nocc + ["e"] * nvir
                sym[i] = "b"
                sym[nocc + a] = "a"
                det = "".join(sym)
                if det in det_to_index:
                    # Use alpha-excited coeff, which is -X_ia (per run_cis)
                    X[i, a] = -coefs[s, det_to_index[det]]
        rho_ex = cis_state_density(mo_coeff, mo_occ, X, state="total")
        td = cis_transition_density(mo_coeff, mo_occ, X)
        exc_energies = cis_result.get("excitation_energies", [0.0] * nstate)
        e_ex = exc_energies[s] if s < len(exc_energies) else 0.0
        states.append(FragmentState(
            label=label, state_index=s,
            energy=mf_energy + e_ex,
            density_total=rho_ex,
            density_alpha=0.5 * rho_ex,
            density_beta=0.5 * rho_ex,
            transition_density=td,
        ))

    return Fragment(
        label=label,
        mol=mol,
        mo_coeff=mo_coeff,
        mo_energy=cis_result.get("mo_energy", np.zeros(nmo)),
        mo_occ=mo_occ,
        states=states,
    )
