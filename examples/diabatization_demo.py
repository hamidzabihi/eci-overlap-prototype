#!/usr/bin/env python3
"""
State-character tracking along a nuclear displacement path.

This demo simulates a photochemically realistic scenario:

  * three electronic states (S0, S1, S2) exist at each geometry
  * along the path, S1 and S2 undergo a sequence of avoided crossings
    with an intervening trivial crossing between S2 and S3
  * at each step, we have an overlap matrix S(R0, R_k) between the
    reference geometry and the current geometry
  * we assign each adiabatic state at R_k to a reference state using
    maximum-overlap matching, and track the resulting "character"

The demo uses SYNTHETIC overlap matrices rather than calling
wfoverlap.x. This is deliberate:

  * The overlap matrices are constructed analytically so the demo
    is fully deterministic and reproducible without any external
    dependency.
  * The state-tracking algorithm is what's being demonstrated here.
    The integration of the algorithm with the wfoverlap engine is
    exercised separately in examples/build_fragment_overlap.py.

See src/diabatization.py for the underlying algorithm.

Usage:
    python examples/diabatization_demo.py [--nsteps 15] [--out plot.png]
"""

from __future__ import annotations
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np

from src.diabatization import (
    track_states, track_reference_state, assignment_table,
    character_matrix, plot_character,
)


# ---------------------------------------------------------------- synthetic path

def rotation_matrix(theta: float) -> np.ndarray:
    """2x2 rotation matrix."""
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s], [s, c]])


def synthetic_overlap(step: int, nsteps: int, nstate: int = 4) -> np.ndarray:
    """
    Build a synthetic (nstate x nstate) overlap matrix that simulates
    the effect of a rotating pair of states along a path.

    The scenario:
      - S0 (ground state) stays pure: S[0,0] ~ 1
      - S1 and S2 rotate through each other: their 2x2 block is a
        rotation by an angle that goes from 0 to pi/2 over the path
      - S3 is a spectator state with moderate overlap to S2

    This produces a path where:
      - at step 0 (reference), S = I
      - at step N-1, S1 and S2 have swapped character
    """
    # Ground state: pure, slight decay
    t = step / max(1, nsteps - 1)
    S = np.zeros((nstate, nstate), dtype=complex)
    S[0, 0] = 1.0 - 0.05 * t   # ~1 throughout

    # S1/S2 rotation: angle from 0 to pi/2
    theta = (np.pi / 2) * t
    R = rotation_matrix(theta)
    S[1:3, 1:3] = R

    # S3: spectator with a small coupling to S2
    S[3, 3] = 1.0 - 0.02 * t
    S[2, 3] = 0.1 * np.sin(np.pi * t)
    S[3, 2] = -S[2, 3]  # antisymmetric

    return S


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nsteps", type=int, default=11,
                    help="number of steps along the path (including reference)")
    ap.add_argument("--nstate", type=int, default=4,
                    help="number of electronic states")
    ap.add_argument("--out", default="diabatization_character.png",
                    help="output plot filename")
    ap.add_argument("--no-plot", action="store_true",
                    help="skip generating the plot")
    args = ap.parse_args()

    nsteps = args.nsteps
    nstate = args.nstate

    print("=" * 64)
    print("Diabatization / state tracking demo (synthetic overlaps)")
    print("=" * 64)
    print(f"  states: {nstate}")
    print(f"  steps:  {nsteps}")
    print()
    print("Scenario: S1 and S2 rotate through each other along the path;")
    print("          S3 is a spectator with a small coupling to S2.")
    print()

    # Build the synthetic overlap matrices
    S_list = [synthetic_overlap(k, nsteps, nstate) for k in range(nsteps)]

    # Print the raw overlap at a few key steps
    print("Reference overlap matrices at key steps:")
    for k in [0, nsteps // 2, nsteps - 1]:
        print(f"\n--- step {k} ---")
        np.set_printoptions(precision=3, suppress=True, linewidth=100)
        print(S_list[k].real)

    # Run the assignment
    assignments = track_states(S_list, threshold=0.3)

    print()
    print("=" * 64)
    print("Assignment table per step (new state -> reference state)")
    print("=" * 64)
    for k, a in enumerate(assignments):
        t = k / max(1, nsteps - 1)
        print(f"\nStep {k:>2}  (t = {t:.3f}):")
        print(assignment_table(a))

    # Character map for each reference state
    print()
    print("=" * 64)
    print("Character maps (|S_ref_i, S_new_j|^2)")
    print("=" * 64)
    for ref in range(nstate):
        char = character_matrix(S_list, reference_index=ref)
        print(f"\nReference state S{ref + 1}:")
        for k in range(nsteps):
            row = " ".join(f"{char[k, j]:6.3f}" for j in range(nstate))
            print(f"  step {k:>2}  [ {row} ]")

    # Plot
    if not args.no_plot:
        try:
            char = character_matrix(S_list, reference_index=1)  # S1, the interesting one
            path_coords = np.linspace(0, 1, nsteps)
            plot_character(
                char,
                path_coords=path_coords,
                reference_label="S2 (reference)",
                outpath=args.out,
            )
            print(f"\nSaved character plot of S2 to {args.out}")
        except ImportError:
            print("\nmatplotlib not installed, skipping plot.")

    # Final verification: at the last step, the assignment should have
    # swapped S1 and S2 (due to the 90-degree rotation)
    a_last = assignments[-1]
    print()
    print("=" * 64)
    print("Verification: state ordering at the end of the path")
    print("=" * 64)
    swapped = (a_last.mapping.get(1) == 2 and a_last.mapping.get(2) == 1)
    if swapped:
        print("  PASS: S1 and S2 have swapped character, as expected")
        print("        from the pi/2 rotation.")
    else:
        print(f"  INFO: assignment at last step: {a_last.mapping}")


if __name__ == "__main__":
    main()
