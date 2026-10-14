"""
Test GS-LE coupling symmetry for the symmetric ethylene dimer.

For a symmetric dimer (A and B identical, at equal distances from
the dimer center), the GS-LE coupling must satisfy:

    H[0,1] = -H[0,2]

The sign is fixed by the gauge choice for the transition densities
(both td vectors have their z-dipole sign fixed in fragment_from_cis_result),
so the magnitudes must match exactly.
"""
import numpy as np
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from examples.build_fragment_overlap import run_cis
from src.excitonic_hamiltonian import (
    build_two_fragment_eci, fragment_from_cis_result,
)
from examples.two_fragment_eci_demo import ethylene_at


def test_gs_le_symmetry():
    sep = 4.0
    geom_A = ethylene_at(z_offset=0.0)
    geom_B = ethylene_at(z_offset=sep)

    r_A = run_cis(geom_A, basis="6-31G", nstates=2)
    r_B = run_cis(geom_B, basis="6-31G", nstates=2)

    frag_A = fragment_from_cis_result("A", r_A)
    frag_B = fragment_from_cis_result("B", r_B)

    H, _ = build_two_fragment_eci(
        frag_A, frag_B,
        include_gs_le=True,
        include_triplets=False,
    )

    V_0A = H[0, 1]
    V_0B = H[0, 2]

    # Symmetric dimer: magnitudes must match to numerical precision
    assert abs(abs(V_0A) - abs(V_0B)) < 1e-8, (
        f"Symmetry broken: |V_0A| = {abs(V_0A):.6e}, "
        f"|V_0B| = {abs(V_0B):.6e}"
    )

    # Physical magnitude: for π-stacked ethylene at 4 Å, expect
    # roughly 1–100 meV
    V_meV = abs(V_0A) * 27211.4
    assert 1.0 < V_meV < 100.0, (
        f"GS-LE magnitude {V_meV:.3f} meV outside expected range"
    )

    # Sanity check on signs
    assert V_0A * V_0B < 0, "Expected opposite signs for symmetric dimer"


def test_gs_le_zero_without_flag():
    """With --include-gs-le=False, the couplings must be exactly zero."""
    sep = 4.0
    geom_A = ethylene_at(z_offset=0.0)
    geom_B = ethylene_at(z_offset=sep)

    r_A = run_cis(geom_A, basis="6-31G", nstates=2)
    r_B = run_cis(geom_B, basis="6-31G", nstates=2)

    frag_A = fragment_from_cis_result("A", r_A)
    frag_B = fragment_from_cis_result("B", r_B)

    H, _ = build_two_fragment_eci(
        frag_A, frag_B,
        include_gs_le=False,
        include_triplets=False,
    )

    assert abs(H[0, 1]) < 1e-12
    assert abs(H[0, 2]) < 1e-12