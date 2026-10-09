# ECI prototype — implementation notes

## Status

- **Singlet ECI**: validated against direct CIS.  Mean absolute
  deviation 122 meV for S1/S2 of the ethylene dimer at 4 Å
  (6-31G basis, TDA).
- **Triplet ECI**: implemented and internally consistent
  (`T2 - T1 = 2|V_TT|` to within 1 meV), but not yet compared
  to a proper diabatic reference.  See "Triplet caveat" below.

## Density construction

The fragment state densities are built in `cis_state_density`.
The correct MO-basis form for a CIS state with amplitudes `X_ia`
is:

    Occ diagonal:  1 - w * sum_a X_ia^2
    Vir diagonal:  w * sum_i X_ia^2
    Occ-vir block: w * X_ia

with `w = 0.5` for a singlet (excitation split equally between
alpha and beta) and `w = 1.0` for the alpha channel of an
M_S = +1 triplet (all excitation in alpha; the beta channel gets
`w = 0`).

The **diagonal occupation-change terms** are essential.  An
earlier implementation only included the off-diagonal coherence
(`cross + cross.T`), which has zero trace and therefore gives
`N_alpha = N_beta` for any state — including a genuine triplet.
That bug is fixed.

## K-term sign rule

The exchange integral K enters the GFT formula with different
signs depending on the physical context:

    Diagonal site energies (state|state):   Y = J - K
    Singlet Frenkel coupling (trans|trans): Y = J + K
    Triplet Frenkel coupling (trans|trans): Y = J - K

The rule is implemented in `Y_integral` via the `spin` keyword
and the `is_transition_F/G` flags.  Getting this wrong shifts
the singlet MAD by ~20 meV and produces unphysical triplet
couplings (`V_TT > V_SS`).

## Triplet caveat

At the RHF/6-31G level, the fragment T1 state of ethylene has
mixed pi->pi* and Rydberg character.  Its CIS amplitude has
`<X_singlet, X_triplet> = 0.475`, i.e. the fragment triplet is
not the same single excitation as the fragment singlet.

Consequence: the ECI `V_TT = +465 meV` is the diabatic Frenkel
coupling between two **mixed-character** fragment triplets.
The direct CIS T1-T2 splitting of 49 meV is the adiabatic
splitting of two **combined-system** triplet eigenstates.  These
are different physical quantities and should not be expected to
match without:

  1. Diabatization of the combined-system CIS triplets onto
     fragment-localized reference states.
  2. Projection of the fragment triplet onto pure pi->pi*
     character (removing the Rydberg admixture).

Both are focused follow-up tasks.

## Files

- `src/excitonic_hamiltonian.py` — ECI Hamiltonian builder,
  fragment density construction, JK integrals, triplet block.
- `examples/build_fragment_overlap.py` — `run_cis` wrapper
  (supports `singlet=True/False`).
- `examples/two_fragment_eci_demo.py` — end-to-end demo with
  direct CIS reference comparison.

## References

- Pitesa, Polonius, Gonzalez, Mai, JCTC 2024, 20, 5609 (ECI)
- Pitesa, Mai, Gonzalez, JPCL 2025, 16, 2800 (RI approx.)
- McWeeny, Methods of Molecular Quantum Mechanics, 2nd ed.

## Sign gauge

Fragment transition densities are put in a fixed phase convention
by enforcing

    mu_z = Tr(td @ Z_int) >= 0

where `Z_int[i,j] = <phi_i | z | phi_j>` is the AO-basis position
operator (PySCF's `int1e_r[2]`).  This is the standard
spectroscopic convention: the transition dipole along the
molecular axis defines the phase.

An earlier implementation used a "largest-magnitude element is
positive" heuristic, which was fragile for multi-configurational
states (mixed pi->pi*/Rydberg) and flipped the sign of both V_SS
and V_TT between runs.  The transition-dipole convention is
reproducible (verified over 5 consecutive runs).

## Test suite

`examples/test_density_spin.py` validates the spin-resolved
fragment densities:

- Singlet:  N_alpha - N_beta = 0
- Triplet:  N_alpha - N_beta = 2   (for M_S = +1)

Both are computed as `Tr(S @ P)` with the AO overlap `S`, which
is the correct electron-count diagnostic in a non-orthogonal
basis.
