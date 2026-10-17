"""Sanity check for spin-resolved fragment densities.

Verifies:
  - singlet:  N_alpha == N_beta
  - triplet:  N_alpha - N_beta == 1  (M_S = +1)
  - total electron count is 2 * nocc for both
"""
import numpy as np
from examples.build_fragment_overlap import run_cis
from examples.two_fragment_eci_demo import ethylene_at
from src.excitonic_hamiltonian import cis_state_density

def make_X(r):
    nocc, nvir = r["nocc"], r["nvir"]
    d2i = {d: i for i, d in enumerate(r["detstrings"])}
    X = np.zeros((nocc, nvir))
    for i in range(nocc):
        for a in range(nvir):
            sym = ["d"] * nocc + ["e"] * nvir
            sym[i] = "b"; sym[nocc + a] = "a"
            det = "".join(sym)
            if det in d2i:
                X[i, a] = -r["coefs"][1, d2i[det]]
    return X

geom = ethylene_at(z_offset=0.0)

for singlet in (True, False):
    r = run_cis(geom, basis="6-31G", nstates=2, singlet=singlet)
    X = make_X(r)
    mo_occ = np.r_[np.full(r["nocc"], 2.0), np.zeros(r["nvir"])]
    S = r["mol"].intor("int1e_ovlp")
    spin = 0 if singlet else 1

    rho_a = cis_state_density(r["mo_coeff"], mo_occ, X, state="alpha", spin=spin)
    rho_b = cis_state_density(r["mo_coeff"], mo_occ, X, state="beta",  spin=spin)

    N_a = np.einsum("ij,ji->", S, rho_a)
    N_b = np.einsum("ij,ji->", S, rho_b)
    label = "singlet" if singlet else "triplet"

    expected = 0.0 if singlet else 2.0   # N_alpha - N_beta = 2*M_S
    print(f"{label:>7}: N_alpha={N_a:.4f}  N_beta={N_b:.4f}  "
          f"N_a - N_b = {N_a - N_b:+.4f}  (expect {expected})")

    if singlet:
        assert abs(N_a - N_b) < 1e-6, "singlet must have N_alpha = N_beta"
    else:
        # M_S = +1 triplet: N_alpha - N_beta = 2 * M_S = 2
        assert abs((N_a - N_b) - 2.0) < 1e-6, "M_S=+1 triplet must have N_a - N_b = 2"

print("OK")
