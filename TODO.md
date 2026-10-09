# Follow-ups

## 1. Rydberg projection on the fragment triplet

Goal: make the fragment T1 a pure pi->pi* state so its ECI
coupling can be compared directly to the combined-system
direct CIS T1-T2 splitting.

Approach:
  a. From `run_cis(..., singlet=False)`, extract the triplet
     amplitudes X_ia for the lowest triplet.
  b. Identify the pi (HOMO) and pi* (LUMO) orbitals from the
     fragment SCF -- e.g. by inspecting the MO coefficients
     near the C=C midpoint, or by the largest contribution to
     the frontier orbitals.
  c. Project X onto the (pi, pi*) subspace and renormalize.
  d. Rebuild the triplet fragment with the projected X and
     recompute V_TT.
  e. Compare to direct CIS T1-T2 splitting.

Success criterion: |V_TT_projected - (delta_TT_CIS / 2)| small,
i.e. the diabatic and adiabatic descriptions converge.

## 2. Diabatization of combined-system triplets

Goal: get the true diabatic coupling V_TT without relying on
fragment-localized states.

Approach:
  a. Compute combined-system CIS triplets (nstates >= 4).
  b. Localize the pi/pi* orbitals on the combined system
     (e.g. Pipek-Mezey localization of the four frontier
     orbitals, or symmetric/antisymmetric combinations of the
     fragment orbitals).
  c. Compute the diabatic states as single-determinant
     configurations in the localized basis:
       T_A = |...; pi_A^up, pi_A*^up>
       T_B = |...; pi_B^up, pi_B*^up>
  d. Transform H from adiabatic to diabatic basis via the
     fragment-reference overlap.
  e. Read V_TT = <T_A|H|T_B> from the transformed matrix.

Success criterion: V_TT_diabatic matches V_TT_ECI to within
the basis-set error.

## 3. Beyond ethylene: a real chromophore

The current demo uses ethylene at 4 A -- a test system.
For the ECI2GAME use case, run the same machinery on a
realistic chromophore dimer (e.g. perylene diimide, or a
small acene) with a proper basis set.  This is where the
prototype actually needs to hold up.

## 4. Performance

The current implementation builds the full cross-fragment
ERI tensor in the AO basis.  For the two-fragment system with
26 AOs each, this is fine.  For larger fragments it will need
the RI approximation from the 2025 paper (JPCL 16, 2800).
