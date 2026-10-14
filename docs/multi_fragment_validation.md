# Multi-Fragment ECI Validation

## What was added

`src/eci_general.py` generalizes the two-fragment ECI construction
to N fragments. Basis structure and math are unchanged; the two-
fragment case is recovered when N = 2.

`examples/n_fragment_test.py` builds N ethylene fragments in a linear
chain, computes the ECI Hamiltonian at FEM level (no GS-LE coupling,
no embedding), and compares the low-lying eigenvalues to direct CIS
on the combined system.

## Results

| N fragments | MAD (meV) | Notes |
|---|---|---|
| 3 | 17.6 | S1-S3 comparison |
| 4 | 19.9 | S1-S4 comparison |

For both tests, the ECI eigenvalues reproduce the direct-CIS
singlet manifold to under 20 meV. This is well below the paper's
FEM-level MAD of 172 meV for G4 (Table 2 of JCTC 2024), because
ethylene is a much simpler system with no strong correlation or
embedding effects.

## Method

- **Basis:** 6-31G
- **Fragments:** ethylene (C2H4), separated by 5 Å along the C=C axis
- **Fragment states:** GS + S1 per fragment (2 states total)
- **Combined-system reference:** TDA-CIS on the full chain, same basis
- **ECI level:** FEM (only GS products and LE products in the basis)

## Correctness

The N-fragment ECI code correctly reproduces the 2-fragment result
when N = 2 (cross-checked against `examples/two_fragment_eci_demo.py`).
The cross-fragment ERIs are built pair-by-pair, and the off-diagonal
Frenkel couplings use the same `_frenkel_coupling` routine as the
two-fragment case.

## What is not validated

- **GS-LE couplings at multi-fragment level.** The multi-fragment
  GS-LE routine (`_gs_le_coupling`) has a TODO for the case where the
  excited fragment's index is greater than the reference fragment.
  The two-fragment case is validated; the N-fragment case is not.
- **Embedding.** The EHF self-consistent ESP cycle is not implemented
  in this prototype. The paper's ECIS + EHF MADs (Table 2, third
  column) are not reproduced.
- **Table 2 G4 reproduction.** Requires the paper's SI geometry
  (def2-SVP basis on the S0-optimized G4 structure). Not attempted
  here; the SI zip is not available locally.

## Next steps

1. Extend `_gs_le_coupling` to handle F > G (the B-first ERI path)
2. Implement the EHF embedding cycle
3. Reproduce Table 2 with the paper's exact geometry and basis
