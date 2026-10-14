#!/usr/bin/env python3
"""
N-fragment ECI validation against direct CIS.

Builds N ethylene fragments in a linear arrangement, runs CIS on each
independently, and compares the ECI excitation energies to direct CIS
on the combined system.

Validation criterion: MAD between ECI S1-S3 and direct CIS S1-S3 should
be under ~150 meV (this is FEM level, no embedding).

Usage:
    python examples/three_fragment_test.py [--sep 4.0]
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
from pyscf import gto, scf, tdscf

from examples.build_fragment_overlap import run_cis
from examples.two_fragment_eci_demo import ethylene_at
from src.excitonic_hamiltonian import fragment_from_cis_result
from src.eci_general import build_eci_general


def direct_cis_energies(atom, basis, nstates):
    mol = gto.M(atom=atom, basis=basis, verbose=0)
    mf = scf.RHF(mol).run()
    td = tdscf.TDA(mf)
    td.nstates = nstates
    td.singlet = True
    td.run()
    return float(mf.e_tot), [float(e) for e in td.e]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sep", type=float, default=4.0,
                    help="fragment-fragment separation (Angstrom)")
    ap.add_argument("--nfrag", type=int, default=3,
                    help="number of fragments in the chain (default 3)")
    ap.add_argument("--basis", default="6-31G")
    ap.add_argument("--nstates", type=int, default=2,
                    help="states per fragment (GS + nstates-1)")
    args = ap.parse_args()

    print("=" * 66)
    print("N-fragment ECI validation")
    print("=" * 66)
    print(f"  separation: {args.sep} A")
    print(f"  basis:      {args.basis}")
    print(f"  nstates/frag: {args.nstates}")
    print()

    # N ethylene fragments in a line
    geoms = [
        ethylene_at(z_offset=i * args.sep) for i in range(args.nfrag)
    ]

    print("Running CIS on each fragment ...")
    frags = []
    for i, g in enumerate(geoms):
        r = run_cis(g, basis=args.basis, nstates=args.nstates)
        frags.append(fragment_from_cis_result(f"F{i}", r))
        print(f"  fragment {i}: SCF = {r['scf_energy']:.6f} Ha, "
              f"excited = {[f'{e*27.2114:.3f}' for e in r['excitation_energies'][1:]]} eV")
    print()

    print("Building ECI Hamiltonian (FEM level) ...")
    H, labels = build_eci_general(frags, include_gs_le=False, verbose=True)
    print()

    eigvals = np.linalg.eigvalsh(H)
    e_gs = eigvals[0]
    eci_ex = [(e - e_gs) * 27.2114 for e in eigvals[1:]]
    eci_sorted = sorted(eci_ex)

    print(f"ECI excitation energies (eV):")
    for i, e in enumerate(eci_sorted[:6]):
        print(f"  S{i+1}: {e:.4f}")

    print()
    print("Running direct CIS on combined system ...")
    combined = []
    for g in geoms:
        combined += g
    e_ref_gs, e_ref_ex = direct_cis_energies(
        combined, args.basis, nstates=args.nstates * args.nfrag
    )
    ref_sorted = sorted([e * 27.2114 for e in e_ref_ex])
    print(f"Direct CIS excitation energies (eV):")
    for i, e in enumerate(ref_sorted[:6]):
        print(f"  S{i+1}: {e:.4f}")

    print()
    print("=" * 66)
    print("Comparison")
    print("=" * 66)
    n = min(6, len(eci_sorted), len(ref_sorted))
    devs = []
    print(f"  {'State':<8}{'ECI (eV)':<14}{'Direct (eV)':<14}{'Dev (meV)':<12}")
    print("  " + "-" * 48)
    for i in range(n):
        dev = (eci_sorted[i] - ref_sorted[i]) * 1000
        devs.append(dev)
        print(f"  S{i+1:<7}{eci_sorted[i]:<14.4f}{ref_sorted[i]:<14.4f}{dev:+.2f}")
    print(f"  MAD: {np.mean(np.abs(devs)):.2f} meV")
    print()

    mad = np.mean(np.abs(devs))
    if mad < 100:
        verdict = "excellent"
    elif mad < 200:
        verdict = "good (FEM level, no embedding)"
    elif mad < 400:
        verdict = "acceptable (large basis mismatch or ordering issue)"
    else:
        verdict = "poor -- check state ordering"
    print(f"  Verdict: {verdict}")


if __name__ == "__main__":
    main()
