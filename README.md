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
| FEM-level ECI (singlet block) | ✅ matches direct CIS (~140 meV MAD) |
| Triplet ECI extension | ✅ validated (T1 dev 78 meV, V_TT ~ 1/R³) |
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
| Triplet Frenkel coupling `V_TT` | 455.81 meV |
| α-only exchange `K` | 17.98 meV |
| T1 (ECI) vs T1 (direct CIS) | 3.6387 eV vs 3.5605 eV (dev +78 meV) |
| Singlet MAD (S1, S2) | 140 meV |
| `2\|V_TT\| = T2 − T1` consistency | holds to < 1 meV |

Separation scan (4–10 Å): `|V_TT|` decays as 1/R³, matching
dipole-dipole theory.

| sep (Å) | V_TT (meV) |
|---|---|
| 4.0 | +455.81 |
| 5.0 | −244.05 |
| 6.0 | +142.16 |
| 8.0 | +60.37 |
| 10.0 | +31.01 |

**Important distinction:** the ECI `V_TT` is a **diabatic** Frenkel
coupling between localized fragment triplet excitations. The direct
CIS T1–T2 splitting is an **adiabatic** splitting of delocalized
combined-system states. These are different observables and agree
only after diabatization of the direct CIS states (see
`NEXT_MILESTONE.md`).

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
