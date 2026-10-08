#!/usr/bin/env python3
"""
Build a fragment-overlap test case from scratch using PySCF.

Computes CIS wavefunctions for H2O at two geometries, converts them
into the three files wfoverlap reads (mos, dets, AO_overl.mixed),
runs wfoverlap.x, and prints the resulting state-overlap matrix.

Usage:
    python examples/build_fragment_overlap.py [--basis 6-31G]

Requires:
    - pyscf
    - wfoverlap.x at ~/sharc4/bin/wfoverlap.x
"""

from __future__ import annotations
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
from pyscf import gto, scf, ci

from src.fragment_overlap import (
    write_mo_file,
    write_dets_file_ci,
    write_ao_overlap_file,
    compute_fragment_overlap,
)


# ---------------------------------------------------------------- geometries

def water_geometry_1():
    """Equilibrium geometry, in Angstrom."""
    return [
        ["O", ( 0.000000,  0.000000,  0.000000)],
        ["H", ( 0.757000,  0.586000,  0.000000)],
        ["H", (-0.757000,  0.586000,  0.000000)],
    ]


def water_geometry_2(displacement=0.05):
    """Slightly displaced geometry (stretch O-H bonds)."""
    g = water_geometry_1()
    scale = 1.0 + displacement
    return [
        [g[0][0], g[0][1]],
        [g[1][0], tuple(x * scale for x in g[1][1])],
        [g[2][0], tuple(x * scale for x in g[2][1])],
    ]


# ---------------------------------------------------------------- PySCF CIS

def run_cis(atom, basis="6-31G", nstates=3, verbose=False):
    """Run RHF + TDA/CIS on the given geometry.

    Returns a dict with the MO coefficients, a common list of detstrings,
    and the CI coefficient matrix for the requested states. The ground
    state is state index 0 (detstring = 'd'*nocc + 'e'*nvir); excited
    singlet states are expanded from the TDA amplitudes X_{ia}.

    The determinant set is the union of {ground state} and, for each
    requested excited state, the singly-excited determinants with
    amplitude > 1e-10.
    """
    from pyscf import tdscf

    mol = gto.M(atom=atom, basis=basis, verbose=0 if not verbose else 4)
    mf = scf.RHF(mol).run()

    nocc = mol.nelectron // 2
    nmo = mf.mo_coeff.shape[1]
    nvir = nmo - nocc
    nao = mf.mo_coeff.shape[0]

    # TDA: excitation amplitudes X[istate][i, a]
    td = tdscf.TDA(mf)
    td.nstates = max(1, nstates - 1)
    td.run()

    # PySCF's td.xy structure (as of the version in use here) is:
    #     td.xy = [ (X_0, Y_0), (X_1, Y_1), ... ]
    # where X_i has shape (nocc, nvir), and Y_i is a placeholder (0 for TDA).
    n_excited = len(td.xy)

    gs_det = "d" * nocc + "e" * nvir

    # Build the coefficient matrix. Column order = sorted list of detstrings.
    # First pass: collect all detstrings and coefficient values.
    det_to_coefs: dict = {}
    det_to_coefs[gs_det] = [1.0] + [0.0] * n_excited

    for istate, (x, y) in enumerate(td.xy):
        # x has shape (nocc, nvir)
        for i in range(nocc):
            for a in range(nvir):
                x_ia = x[i, a]
                if abs(x_ia) < 1e-10:
                    continue
                # Alpha excitation (SHARC convention):
                # the occupied orbital i loses its alpha electron (d -> b),
                # the virtual orbital a gains an alpha electron (e -> a).
                sym_a = ["d"] * nocc + ["e"] * nvir
                sym_a[i] = "b"
                sym_a[nocc + a] = "a"
                det_a = "".join(sym_a)

                # Beta excitation (SHARC convention):
                # the occupied orbital i loses its beta electron (d -> a),
                # the virtual orbital a gains a beta electron (e -> b).
                sym_b = ["d"] * nocc + ["e"] * nvir
                sym_b[i] = "a"
                sym_b[nocc + a] = "b"
                det_b = "".join(sym_b)

                if det_a not in det_to_coefs:
                    det_to_coefs[det_a] = [0.0] * (1 + n_excited)
                if det_b not in det_to_coefs:
                    det_to_coefs[det_b] = [0.0] * (1 + n_excited)

                # Singlet combination, following SHARC_GAUSSIAN's sign convention:
                # alpha excitation gets a minus sign, beta excitation a plus sign.
                #
                # ACCUMULATE (not overwrite): the same detstring can arise from
                # different (i, a) pairs in larger CIS spaces (e.g. formaldehyde),
                # and their contributions must be summed. Overwriting loses
                # normalization and produces coefficients with |c| > 1.
                det_to_coefs[det_a][1 + istate] += -x_ia
                det_to_coefs[det_b][1 + istate] +=  x_ia

    detstrings = sorted(det_to_coefs.keys())
    nstate = 1 + n_excited
    coefs = np.zeros((nstate, len(detstrings)))
    for d, det in enumerate(detstrings):
        coefs[:, d] = det_to_coefs[det]

    # Collect excitation energies from PySCF's TDA result (in Hartree)
    excitation_energies = [0.0]  # GS has zero excitation energy
    for istate in range(n_excited):
        e_ex = float(td.e[istate])
        excitation_energies.append(e_ex)

    return {
        "mol": mol,
        "mo_coeff": mf.mo_coeff,
        "mo_energy": mf.mo_energy,
        "mo_occ": mf.mo_occ,
        "scf_energy": float(mf.e_tot),
        "excitation_energies": excitation_energies,
        "detstrings": detstrings,
        "coefs": coefs,
        "nmo": nmo,
        "nocc": nocc,
        "nvir": nvir,
        "nao": nao,
        "nstate": nstate,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--basis", default="6-31G")
    parser.add_argument("--nstates", type=int, default=3)
    parser.add_argument("--wfoverlap", default=str(Path.home() / "sharc4" / "bin" / "wfoverlap.x"))
    parser.add_argument("--workdir", default="/tmp/eci_overlap_test")
    parser.add_argument("--displacement", type=float, default=0.05)
    args = parser.parse_args()

    workdir = Path(args.workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    print(f"Working directory: {workdir}")
    print(f"Basis: {args.basis}  States: {args.nstates}  Displacement: {args.displacement}")

    # Geometry A
    atom_a = water_geometry_1()
    print("\nComputing CIS at geometry A...")
    result_a = run_cis(atom_a, basis=args.basis, nstates=args.nstates)
    print(f"  nAO={result_a['nao']}  nMO={result_a['nmo']}  nstates={args.nstates}")

    # Geometry B
    atom_b = water_geometry_2(displacement=args.displacement)
    print("\nComputing CIS at geometry B...")
    result_b = run_cis(atom_b, basis=args.basis, nstates=args.nstates)
    print(f"  nAO={result_b['nao']}  nMO={result_b['nmo']}  nstates={args.nstates}")

    # Write MO files
    write_mo_file(workdir / "mos_a", result_a["mo_coeff"])
    write_mo_file(workdir / "mos_b", result_b["mo_coeff"])

    # Write dets files
    # Note: the two geometries might have different determinant lists if
    # the CI vectors differ. For validation, we use the union of determinants.
    all_dets = sorted(set(result_a["detstrings"]) | set(result_b["detstrings"]))
    det_index = {d: i for i, d in enumerate(all_dets)}
    nstate_a = result_a["nstate"]
    nstate_b = result_b["nstate"]
    if nstate_a != nstate_b:
        raise RuntimeError(f"states differ: {nstate_a} vs {nstate_b}")
    nstate = nstate_a
    coefs_a = np.zeros((nstate, len(all_dets)))
    coefs_b = np.zeros((nstate, len(all_dets)))
    for i, d in enumerate(result_a["detstrings"]):
        coefs_a[:, det_index[d]] = result_a["coefs"][:, i]
    for i, d in enumerate(result_b["detstrings"]):
        coefs_b[:, det_index[d]] = result_b["coefs"][:, i]

    write_dets_file_ci(workdir / "dets_a", all_dets, coefs_a, result_a["nmo"])
    write_dets_file_ci(workdir / "dets_b", all_dets, coefs_b, result_b["nmo"])

    # Write AO overlap
    mol_old = result_a["mol"]
    mol_new = result_b["mol"]
    mol_conc = gto.conc_mol(mol_old, mol_new)
    mol_conc.build()
    nao_old = mol_old.nao
    sao_mixed = mol_conc.intor("int1e_ovlp")[:nao_old, nao_old:]
    write_ao_overlap_file(workdir / "AO_overl.mixed", sao_mixed)

    print(f"\nWrote input files to {workdir}")
    for f in ["mos_a", "mos_b", "dets_a", "dets_b", "AO_overl.mixed"]:
        print(f"  {f}  ({ (workdir / f).stat().st_size } bytes)")

    # Run wfoverlap
    print("\nRunning wfoverlap.x...")
    S = compute_fragment_overlap(
        mo_a="mos_a",
        mo_b="mos_b",
        det_a="dets_a",
        det_b="dets_b",
        ao_overlap="AO_overl.mixed",
        workdir=workdir,
        wfoverlap_exe=args.wfoverlap,
        keep_inputs=True,
    )

    print(f"\nOverlap matrix S (shape {S.shape}):")
    print(S)
    if S.shape[0] == S.shape[1]:
        print("\n|S| close to identity?", np.allclose(np.abs(S), np.eye(S.shape[0]), atol=0.2))
    else:
        print("\nShape is not square -- skipping identity check.")


if __name__ == "__main__":
    main()
