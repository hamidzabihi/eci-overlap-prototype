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

## Status

| Stage | Status |
|---|---|
| Fragment overlap primitive | ✅ tested |
| ESD overlap product formula | ✅ tested |
| ECSF contraction / state overlaps | ✅ tested |
| Wavefunction phase correction | ✅ tested |
| Diabatization / state tracking | ✅ tested |
| ECI Hamiltonian from fragments | ✅ tested |
| **FEM-level ECI** | ✅ matches direct CIS (~140 meV MAD) |
| **ECIS-level ECI (with GS-LE)** | ❌ FEM level only; GS-LE not yet validated |
| Full SHARC_ECI.py integration | ❌ future work |
| Energy gradients | ❌ future work |
| Spin-orbit couplings | ❌ future work |

The remaining stages map directly onto the ECI2GAME postdoc deliverables.

## Known limitations

See `NOTES.md` for a detailed writeup, including the `wfoverlap.x`
build issue (gfortran + -i4 LAPACK causes silent numerical corruption
for the native-ascii input path).

## Related work

- SHARC 4 source: https://github.com/sharc-md/sharc4
- `cis_nto` (CP2K reader PR): https://github.com/marin-sapunar/cis_nto/pull/4
