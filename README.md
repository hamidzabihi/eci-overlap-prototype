# eci-overlap-prototype

A prototype of the **ECI → wfoverlap → SHARC-NAMD** pipeline.

The repo has two complementary goals:

1. **Overlap pipeline:** compute the state-overlap matrix `S_ij` for
   an ECI wavefunction at two geometries, which is what SHARC's
   surface-hopping NAMD needs at every time step.
2. **ECI method demo:** a fragment-based implementation of the
   Excitonic Configuration Interaction method (Piteša et al., JCTC
   2024; JPCL 2025) that reproduces the excitonic splitting for a
   two-chromophore system.

## Quick start

    pip install numpy pytest pyscf matplotlib
    pytest -v                                       # 88 tests
    python examples/build_fragment_overlap.py       # fragment overlap from PySCF
    python examples/two_fragment_eci_demo.py        # working ECIS calculation
    python examples/diabatization_demo.py           # state tracking (synthetic)
    python examples/gap_map_demo.py                 # QM.out section gap map
    python examples/pipeline_report.py              # pipeline stage status

## Requirements

- Python 3.10+, NumPy, pytest
- PySCF for the demos (not required for the module tests)
- Matplotlib for the diabatization plot (optional)
- A clone of `sharc-md/sharc4` at `~/sharc4` with a compiled
  `wfoverlap.x` (see `BUILDING.md`)

## Layout

    src/
      run_wfoverlap.py         # wrapper around the SHARC wfoverlap binary
      qmout_parser.py          # parser for SHARC's QM.out format + gap map
      esd_to_dets.py           # ESD -> detstring mapper
      fragment_overlap.py      # fragment file writers + wfoverlap invocation
      esd_overlap.py           # ESD overlap product formula
      state_overlap.py         # ECSF contraction + state overlaps
      phase_tracking.py        # wavefunction phase correction (simple + robust)
      diabatization.py         # state tracking across a geometry path
      excitonic_hamiltonian.py # ECI Hamiltonian from fragment calculations
    tests/                     # 88 tests, all passing
    examples/                  # runnable demos
    BUILDING.md                # how to build wfoverlap with gfortran
    NOTES.md                   # technical writeup and known limitations
    NEXT_MILESTONE.md          # proposed diabatization validation milestone

## Status

| Stage | Status |
|---|---|
| Fragment overlap primitive | ✅ tested |
| ESD overlap product formula | ✅ tested |
| ECSF contraction / state overlaps | ✅ tested |
| Wavefunction phase correction | ✅ tested |
| Diabatization / state tracking | ✅ tested |
| ECI Hamiltonian from fragments | ✅ tested |
| FEM-level ECI (singlet block) | ✅ matches direct CIS (122 meV MAD) |
| Triplet ECI extension | ⚠️ implemented; diabatic-vs-adiabatic caveat (see NOTES.md) |
| ECIS-level ECI (with GS-LE) | ❌ FEM level only; GS-LE not yet validated |
| Full SHARC_ECI.py integration | ❌ future work |
| Energy gradients | ❌ future work |
| Spin-orbit couplings | ❌ future work |

The remaining stages map directly onto the ECI2GAME postdoc deliverables.

## Recent milestone — triplet extension (tag `triplet-extension-v1`)

The ECI Hamiltonian now includes the triplet block alongside the
singlet block, with a 6×6 basis:

    [ GS, S_A, S_B | T_A, T_B, T_A-T_B ]

Key results at 4 Å (ethylene dimer, 6-31G):

| Quantity | Value |
|---|---|
| Triplet Frenkel coupling `V_TT` | −464.80 meV |
| α-only exchange `K` | 17.98 meV |
| T1 (ECI) vs T1 (direct CIS) | 3.1236 eV vs 3.5605 eV (dev −437 meV) |
| Singlet MAD (S1, S2) | 122 meV |
| `2\|V_TT\| = T2 − T1` consistency | holds to < 1 meV |

Separation scan (4–10 Å, singlet Frenkel coupling — the validated
quantity):

| sep (Å) | V_SS (meV) |
|---|---|
| 4.0  | −299.71 |
| 5.0  | −158.85 |
| 6.0  |  −96.25 |
| 8.0  |  −42.62 |
| 10.0 |  −22.31 |

Log-log fit over 4–10 Å: `|V_SS| ~ R^(−2.84)`, approaching the
`R^(−3)` asymptote expected for dipole-dipole coupling.

**Important distinction:** the ECI `V_TT` is a **diabatic** Frenkel
coupling between localized fragment triplet excitations. The direct
CIS T1–T2 splitting is an **adiabatic** splitting of delocalized
combined-system states. These are different observables and agree
only after (a) diabatization of the direct CIS states, and (b)
verification that the fragment triplet and the combined-system
triplet share the same orbital character.

At the RHF/6-31G level, the fragment T1 of ethylene has mixed
`pi -> pi*` and Rydberg character (`<X_singlet, X_triplet> = 0.475`),
so (b) does not hold and the ECI `V_TT` should not be expected to
match the direct CIS T1–T2 splitting. See `NOTES.md` for the full
caveat and `TODO.md` for the follow-up plan (Rydberg projection or
diabatization).

The ECI assumes strongly orthogonal fragments (no inter-fragment
density overlap). For ethylene this holds above ~3.5 Å; the demo
warns below that threshold. See `NOTES.md` for the full validation
summary.

## Known limitations

See `NOTES.md` for a detailed writeup, including:

- The `wfoverlap.x` build issue (gfortran + `-i4` LAPACK causes
  silent numerical corruption for the native-ascii input path).
- The GS-LE coupling is not yet validated; the demo runs at the FEM
  level (`--include-gs-le` exists but the code path is documented
  as unreliable).
- The ECI is only valid for non-overlapping fragments
  (separation ≳ 3.5 Å for the ethylene demo).

## Related work

- SHARC 4 source: https://github.com/sharc-md/sharc4
- `cis_nto` (CP2K reader PR): https://github.com/marin-sapunar/cis_nto/pull/4
