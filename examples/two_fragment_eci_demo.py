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
    fragment_from_cis_result,
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

def direct_cis_energies(atom, basis, nstates):
    """Run CIS on the full system and return the low-lying excitation energies."""
    mol = gto.M(atom=atom, basis=basis, verbose=0)
    mf = scf.RHF(mol).run()
    td = tdscf.TDA(mf)
    td.nstates = nstates
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
    r_A = run_cis(geom_A, basis=args.basis, nstates=args.nstates)
    print(f"  scf_energy:  {r_A['scf_energy']:.8f} Ha")
    print(f"  excitations: {[f'{e:.5f}' for e in r_A['excitation_energies']]}")

    print("Running CIS on fragment B ...")
    r_B = run_cis(geom_B, basis=args.basis, nstates=args.nstates)
    print(f"  scf_energy:  {r_B['scf_energy']:.8f} Ha")
    print(f"  excitations: {[f'{e:.5f}' for e in r_B['excitation_energies']]}")
    print()

    # --- Step 2: build Fragment objects
    frag_A = fragment_from_cis_result("A", r_A)
    frag_B = fragment_from_cis_result("B", r_B)

    # --- Step 3: assemble the ECI Hamiltonian
    print("Assembling ECI Hamiltonian ...")
    H, labels = build_two_fragment_eci(
        frag_A, frag_B,
        include_gs_le=args.include_gs_le,
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

    # --- Step 5: direct reference calculation
    geom_full = combined_geometry(args.separation)
    print("Running direct CIS on combined two-ethylene system (reference) ...")
    e_ref_gs, e_ref_ex = direct_cis_energies(geom_full, args.basis, args.nstates * 2)

    print(f"  SCF energy:  {e_ref_gs:.8f} Ha")
    print(f"  Excitations (eV): {[f'{e * 27.2114:.5f}' for e in e_ref_ex]}")
    print()

    # --- Step 6: comparison
    print("=" * 66)
    print("Comparison: ECI vs direct CIS")
    print("=" * 66)
    print(f"{'State':<8} {'ECI (eV)':<14} {'Direct (eV)':<14} {'Deviation (meV)':<16}")
    print("-" * 52)

    n_compare = min(len(eigvals) - 1, len(e_ref_ex))
    deviations = []
    for i in range(n_compare):
        e_eci = (eigvals[i + 1] - e_gs) * 27.2114
        e_dir = e_ref_ex[i] * 27.2114
        dev = (e_eci - e_dir) * 1000  # meV
        deviations.append(dev)
        print(f"S{i + 1:<7} {e_eci:<14.5f} {e_dir:<14.5f} {dev:+.2f}")

    print()
    if deviations:
        mad = np.mean(np.abs(deviations))
        print(f"Mean absolute deviation: {mad:.2f} meV")
    print()

    # --- Step 7: physical interpretation
    print("=" * 66)
    print("Interpretation")
    print("=" * 66)

    # At large separation, the coupling should vanish; at small separation,
    # the LE states should split into symmetric/antisymmetric combinations.
    if len(eigvals) >= 3:
        splitting = (eigvals[2] - eigvals[1]) * 27.2114
        print(f"Excitonic splitting (E_2 - E_1): {splitting * 1000:.2f} meV")
        if abs(splitting) < 0.01:
            print("  -> negligible coupling (fragments effectively isolated)")
        else:
            print("  -> nonzero coupling: ECI captures inter-fragment excitonic interaction")


if __name__ == "__main__":
    main()
