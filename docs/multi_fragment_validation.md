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

MAD against direct CIS, as a function of fragment separation:

| sep (Å) | 3-frag MAD (meV) | 4-frag MAD (meV) | Note |
|---------|-----------------|-----------------|------|
| 4.0 | 195.7 | 220.1 | fragments overlap; strong-orthogonality assumption breaks |
| 5.0 |  17.6 |  19.9 | clean non-overlapping regime |
| 6.0 |   3.1 |   3.5 | asymptotic regime |
| 8.0 |   0.5 |   0.6 | essentially exact |

The absolute accuracy depends strongly on separation. At 5 Å and
beyond, the ECI eigenvalues reproduce the direct-CIS singlet
manifold to under 20 meV — well below the paper's FEM-level MAD
of 172 meV for G4 (Table 2 of JCTC 2024). This is because ethylene
is a much simpler system (no strong correlation, no embedding).
At 4 Å, the ethylene fragments are close enough to overlap
significantly, and the ECI strong-orthogonality assumption fails.

The physical diagnostic is `Y_LE,LE|00`, the two-electron integral
of the excited-state density of one fragment with the ground-state
density of another:

| sep (Å) | Y_LE,LE|00 (Ha) | (meV) |
|---------|----------------|-------|
| 4.0 | -0.04420060 | -1202 |
| 5.0 | -0.00078901 |   -21 |
| 6.0 | +0.00004426 |    +1 |
| 8.0 | +0.00000241 |    +0.07 |

Below ~5 Å this integral is of order 100 meV or more, indicating
that the two fragments are not independent. Above 5 Å it drops
rapidly toward zero, and the ECI error follows.

## Method

- **Basis:** 6-31G
- **Fragments:** ethylene (C2H4), separated by 5 Å along the C=C axis
  (default in `n_fragment_test.py`; use `--sep` to change)
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
