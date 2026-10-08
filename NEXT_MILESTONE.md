# Next milestone: Diabatization of direct CIS for ECI validation

## Motivation

The triplet extension (tag: triplet-extension-v1) computes the
DIABATIC Frenkel coupling between localized fragment excitations:

    V_TT(ECI) = 455.81 meV at 4 Angstrom (ethylene, 6-31G)

The direct CIS T1-T2 splitting at the same geometry is only
48.6 meV -- but this is the ADIABATIC splitting of delocalized
combined-system states, not the diabatic coupling.  The two are
different observables.

To do a genuine apples-to-apples validation, we need to diabatize
the direct CIS states onto the fragment-localized basis and
extract the diabatic coupling matrix element <T_A|H|T_B>.

## Approach

1. Run CIS on the combined two-fragment system (both singlet and
   triplet manifolds).

2. Localize the combined-system MOs onto fragments.  Options:
   - Boys localization (maximizes spatial separation)
   - Pipek-Mezey (maximizes atomic-orbital locality)
   - Fragment-orbital projection (project fragment MOs onto the
     combined-system AO basis)

3. Re-express the TDA eigenvectors in the localized MO basis.

4. Construct the diabatic states as localized excitations:
       |T_A> = local excitation on fragment A
       |T_B> = local excitation on fragment B

5. Extract the diabatic Hamiltonian:
       <T_A|H|T_B> = V_TT(diabatic)

6. Compare to the ECI V_TT across a separation scan.

## Expected outcome

- V_TT(ECI) and V_TT(diabatized CIS) should agree to within the
  method difference (ECI uses fragment reference calculations;
  diabatized CIS uses the combined-system reference).
- Agreement to ~10-20% would validate the ECI triplet coupling.
- The 1/R^3 scaling should appear in both.

## Implementation notes

- PySCF has `pyscf.lo.Boys` and `pyscf.lo.PM` for localization.
- For fragment-orbital projection, the existing `one_electron_operator`
  and `nuclear_terms` machinery in `src/excitonic_hamiltonian.py`
  already handles the mixed fragment basis -- the same idea
  extends to MO localization.
- The `run_cis` function already returns the full CIS vectors;
  re-expressing them in the localized basis is a unitary
  transformation.

## Success criteria

1. V_TT(diabatized CIS) matches V_TT(ECI) to within 20% across
   the 4-10 Angstrom range.
2. Both show 1/R^3 decay.
3. The T1 energies agree to within 100 meV.

If (1) and (2) hold, the ECI triplet extension is fully validated
against an independent diabatic reference.

## Estimated effort

- Localization + re-expression: ~2-3 hours of coding
- Validation scan + analysis: ~1 hour
- Writeup: ~1 hour

A well-scoped next milestone.
