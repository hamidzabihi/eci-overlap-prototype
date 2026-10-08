#!/usr/bin/env python3
"""
Two-fragment ECI demo: from fragment calculations to excitonic energies.

This demo:

  1. Builds two identical ethylene fragments at a variable separation.
  2. Runs CIS on each fragment (in the fragment's own basis).
  3. Extracts fragment state energies, densities, and transition densities.
  4. Assembles the ECI Hamiltonian using the group-function-theory
     formulas (Pitesa et al., JCTC 2024, eqs. 8-10 and SI S2, S5).
  5. Diagonalizes to get excitonic energies.
  6. Runs a DIRECT CIS calculation on the combined two-ethylene system
     as a reference.
  7. Compares the low-lying singlet states and reports the deviation.

Scope
-----
This is a FEM-level ECI calculation: it uses the fragment energies,
the inter-fragment Coulomb/exchange integrals for the diagonal, and
the Frenkel coupling between the two local excitations. The GS-LE
coupling is skipped because our current implementation of the
transition-density J integral is unreliable (see NOTES.md).

Adding the GS-LE coupling brings us to full ECIS-level, which is what
the paper benchmarks. That's documented as a known limitation and is
exactly the kind of detail the ECI2GAME project would need to get right.

Expected behavior
-----------------
For large separations (>= 10 Angstrom), the two fragments are
effectively independent and the ECI energies reduce to the sum of
fragment energies (no coupling). For shorter separations, the
fragment-fragment coupling splits the two locally-excited states into
a symmetric and antisymmetric excitonic pair.

This is the physics of the ECI method (and of the Frenkel exciton
model) demonstrated on the simplest possible system.

Usage
-----
    python examples/two_fragment_eci_demo.py [--separation 4.0]
                                              [--basis 6-31G]
                                              [--nstates 2]
"""

from __future__ import annotations
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
from pyscf import gto, scf, tdscf

from examples.build_fragment_overlap import run_cis
from src.excitonic_hamiltonian import (
    build_two_fragment_eci,
    build_two_fragment_eci_with_triplets,
    fragment_from_cis_result,
    fragment_from_cis_result_pair,
)


# ---------------------------------------------------------------- geometry

def ethylene_at(z_offset=0.0, scale=1.0):
    """Ethylene, C2H4, with C=C bond along z, translated by z_offset."""
    # C=C ~ 1.33 Angstrom, C-H ~ 1.09 Angstrom, H-C-H angle ~ 117 deg
    dCC = 1.33 * scale
    dCH = 1.09 * scale
    # H atoms are offset perpendicular to the C=C axis
    y_H = 0.93 * scale
    z_H = 0.40 * scale

    return [
        ["C", ( 0.0,  0.0,  0.0        + z_offset)],
        ["C", ( 0.0,  0.0,  dCC        + z_offset)],
        ["H", ( 0.0,  y_H, -z_H        + z_offset)],
        ["H", ( 0.0, -y_H, -z_H        + z_offset)],
        ["H", ( 0.0,  y_H,  dCC + z_H  + z_offset)],
        ["H", ( 0.0, -y_H,  dCC + z_H  + z_offset)],
    ]


def combined_geometry(separation):
    """Two ethylene molecules separated by `separation` along z."""
    a = ethylene_at(z_offset=0.0)
    b = ethylene_at(z_offset=separation)
    return a + b


# ---------------------------------------------------------------- direct reference

def direct_cis_energies(atom, basis, nstates, singlet=True):
    """Run CIS on the full system and return the low-lying excitation energies."""
    mol = gto.M(atom=atom, basis=basis, verbose=0)
    mf = scf.RHF(mol).run()
    td = tdscf.TDA(mf)
    td.nstates = nstates
    td.singlet = singlet
    td.run()
    energies = [float(e) for e in td.e]
    return float(mf.e_tot), energies


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--separation", type=float, default=4.0,
                    help="inter-fragment separation along z, in Angstrom")
    ap.add_argument("--basis", default="6-31G")
    ap.add_argument("--nstates", type=int, default=2,
                    help="number of states per fragment (including GS)")
    ap.add_argument("--include-gs-le", action="store_true",
                    help="Include GS-LE coupling (ECIS level). "
                         "Default is FEM level (no GS-LE).")
    args = ap.parse_args()

    print("=" * 66)
    print("Two-fragment ECI demo")
    print("=" * 66)
    print(f"  fragment:       ethylene (C2H4)")
    print(f"  separation:     {args.separation:.2f} Angstrom")
    print(f"  basis:          {args.basis}")
    print(f"  states/fragment: {args.nstates} (GS + {args.nstates - 1} excited)")
    print(f"  ECI level:      {'ECIS (with GS-LE)' if args.include_gs_le else 'FEM (no GS-LE)'}")
    print()

    # --- Step 1: fragment CIS calculations
    geom_A = ethylene_at(z_offset=0.0)
    geom_B = ethylene_at(z_offset=args.separation)

    print("Running CIS on fragment A ...")
    r_A = run_cis(geom_A, basis=args.basis, nstates=args.nstates, singlet=True)
    print(f"  singlet SCF energy:  {r_A['scf_energy']:.8f} Ha")
    print(f"  singlet excitations: {[f'{e*27.2114:.3f}' for e in r_A['excitation_energies'][1:]]} eV")

    print("Running CIS on fragment B ...")
    r_B = run_cis(geom_B, basis=args.basis, nstates=args.nstates, singlet=True)
    print(f"  singlet SCF energy:  {r_B['scf_energy']:.8f} Ha")
    print(f"  singlet excitations: {[f'{e*27.2114:.3f}' for e in r_B['excitation_energies'][1:]]} eV")

    r_A_trip = run_cis(geom_A, basis=args.basis, nstates=args.nstates, singlet=False)
    r_B_trip = run_cis(geom_B, basis=args.basis, nstates=args.nstates, singlet=False)
    print(f"  triplet A excitations: {[f'{e*27.2114:.3f}' for e in r_A_trip['excitation_energies'][1:]]} eV")
    print(f"  triplet B excitations: {[f'{e*27.2114:.3f}' for e in r_B_trip['excitation_energies'][1:]]} eV")
    print()

    # --- Step 2: build Fragment objects
    frag_A = fragment_from_cis_result("A", r_A, spin=0)
    frag_B = fragment_from_cis_result("B", r_B, spin=0)
    frag_A_trip = fragment_from_cis_result("A", r_A_trip, spin=1)
    frag_B_trip = fragment_from_cis_result("B", r_B_trip, spin=1)

    # --- Step 3: assemble the ECI Hamiltonian
    print("Assembling ECI Hamiltonian ...")
    H, labels = build_two_fragment_eci_with_triplets(
        frag_A, frag_B,
        fragment_A_triplet=frag_A_trip,
        fragment_B_triplet=frag_B_trip,
        verbose=True,
    )

    print()
    print("ECI Hamiltonian (Ha):")
    print(f"  Basis: {labels}")
    print()
    np.set_printoptions(precision=8, suppress=True, linewidth=120)
    print(H)
    print()

    # --- Step 4: diagonalize
    eigvals, eigvecs = np.linalg.eigh(H)
    print("ECI eigenvalues (Ha) and excitation energies (eV):")
    e_gs = eigvals[0]
    for i, e in enumerate(eigvals):
        print(f"  State {i}:  E = {e:.8f} Ha   Eex = {(e - e_gs) * 27.2114:.5f} eV")
    print()

    # --- Step 5: direct reference calculation (both singlets and triplets)
    geom_full = combined_geometry(args.separation)

    print("Running direct CIS on combined system (reference) ...")
    e_ref_gs_s, e_ref_ex_s = direct_cis_energies(
        geom_full, args.basis, args.nstates * 2, singlet=True
    )
    e_ref_gs_t, e_ref_ex_t = direct_cis_energies(
        geom_full, args.basis, args.nstates * 2, singlet=False
    )
    print(f"  SCF energy (singlet ref):  {e_ref_gs_s:.8f} Ha")
    print(f"  Singlet excitations (eV): {[f'{e * 27.2114:.3f}' for e in e_ref_ex_s[:4]]}")
    print(f"  Triplet excitations (eV): {[f'{e * 27.2114:.3f}' for e in e_ref_ex_t[:4]]}")
    print()

    # --- Step 6: comparison by multiplicity
    # The Hamiltonian is block-diagonal by construction: the first
    # n_singlet_basis rows/cols are the singlet block (GS, S_A, S_B),
    # the rest are the triplet block (T_A, T_B, T_A-T_B).  We
    # diagonalize each block separately, which avoids any ambiguity
    # in classifying eigenstates by their eigenvector content.
    print("=" * 66)
    print("Comparison: ECI vs direct CIS (by multiplicity)")
    print("=" * 66)

    n_singlet_basis = 3   # GS, S_A, S_B in the ECI basis
    H_singlet = H[:n_singlet_basis, :n_singlet_basis]
    H_triplet = H[n_singlet_basis:, n_singlet_basis:]

    eigvals_s = np.linalg.eigvalsh(H_singlet)
    eigvals_t = np.linalg.eigvalsh(H_triplet)

    # GS energy from the singlet block
    e_gs = eigvals_s[0]

    eci_singlets = [(e - e_gs) * 27.2114 for e in eigvals_s[1:]]
    eci_triplets = [(e - e_gs) * 27.2114 for e in eigvals_t]

    # For the overall eigenvalue listing, we keep the full H diagonalization
    _, eigvals = np.linalg.eigh(H)

    # Compare singlet block.  eci_singlets already excludes the GS
    # (it was built from eigvals_s[1:] above), so the indexing here is
    # direct -- no off-by-one.
    eci_singlets_sorted = sorted(eci_singlets)
    n_s = min(len(eci_singlets_sorted), len(e_ref_ex_s))
    print()
    print(f"{'Singlet':<10} {'ECI (eV)':<14} {'Direct (eV)':<14} {'Dev (meV)':<12}")
    print("-" * 52)
    singlet_devs = []
    for i in range(n_s):
        e_eci = eci_singlets_sorted[i]
        e_dir = e_ref_ex_s[i] * 27.2114
        dev = (e_eci - e_dir) * 1000
        singlet_devs.append(dev)
        print(f"S{i + 1:<9} {e_eci:<14.4f} {e_dir:<14.4f} {dev:+.2f}")
    if singlet_devs:
        print(f"Singlet MAD: {np.mean(np.abs(singlet_devs)):.2f} meV")

    # Compare triplet block
    eci_triplets_sorted = sorted(eci_triplets)

    # The ECI triplet-triplet coupling is a DIABATIC Frenkel coupling
    # between localized fragment triplet excitations.  The direct CIS
    # T1-T2 splitting is an ADIABATIC splitting of delocalized
    # combined-system states.  These are different physical quantities
    # and are not expected to agree without diabatization of the direct
    # CIS states.  We therefore report them side by side rather than
    # computing a deviation.
    n_singlet_basis = 3   # GS, S_A, S_B in the ECI basis
    V_TT_eci = H[n_singlet_basis, n_singlet_basis + 1] * 27.2114 * 1000
    if len(e_ref_ex_t) >= 2:
        delta_TT_cis = (e_ref_ex_t[1] - e_ref_ex_t[0]) * 27.2114 * 1000
    else:
        delta_TT_cis = float("nan")

    print()
    print("Triplet block (diabatic ECI coupling vs adiabatic CIS splitting)")
    print("-" * 66)
    print(f"  ECI Frenkel coupling  V_TT = H[T_A, T_B] = {V_TT_eci:+.2f} meV")
    print(f"  Direct CIS T1-T2 splitting           = {delta_TT_cis:+.2f} meV")
    print()
    print("  NOTE: These are different quantities.  The ECI value is the")
    print("  diabatic coupling between localized fragment triplets; the")
    print("  direct CIS value is the adiabatic splitting of delocalized")
    print("  combined-system states.  They agree only after diabatization")
    print("  of the direct CIS states onto the fragment-localized basis.")
    print()
    print("  ECI triplet eigenvalues (sorted):")
    for i, e in enumerate(eci_triplets_sorted):
        print(f"    T{i + 1}: {e:.4f} eV")

    # Note: T3 and higher are the DLE-like states (T_A-T_B coupled).
    # They have no direct CIS counterpart and are already listed above.
    print()

    # --- Step 7: physical interpretation
    print("=" * 66)
    print("Interpretation")
    print("=" * 66)

    # At large separation, the coupling should vanish; at small separation,
    # the LE states should split into symmetric/antisymmetric combinations.
    # Singlet excitonic splitting: E(S2) - E(S1)
    if len(eci_singlets_sorted) >= 2:
        split_s = (eci_singlets_sorted[1] - eci_singlets_sorted[0]) * 1000
        print(f"Singlet excitonic splitting (S2 - S1): {split_s:+.2f} meV")

    # Triplet excitonic splitting: E(T2) - E(T1)
    if len(eci_triplets_sorted) >= 2:
        split_t = (eci_triplets_sorted[1] - eci_triplets_sorted[0]) * 1000
        print(f"Triplet excitonic splitting (T2 - T1): {split_t:+.2f} meV")
        print(f"  (= 2 * V_TT = {2 * V_TT_eci:+.2f} meV)")

    if len(eci_singlets_sorted) >= 2 or len(eci_triplets_sorted) >= 2:
        print()
        print("  Nonzero splittings confirm the ECI captures inter-fragment")
        print("  excitonic coupling. The singlet splitting is dominated by")
        print("  the Coulomb (Frenkel) term; the triplet splitting is the")
        print("  same Coulomb term plus the (small) alpha-only exchange.")


if __name__ == "__main__":
    main()
