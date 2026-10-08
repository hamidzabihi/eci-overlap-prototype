# eci-overlap-prototype

A prototype of the **ECI → wfoverlap → SHARC-NAMD** overlap pipeline.

This repo implements the chain of operations needed to compute the
**state-overlap matrix S_ij** for an ECI wavefunction at two geometries,
which is the quantity SHARC's surface-hopping NAMD machinery needs at
every time step.

## What this actually does

Given two geometries of a molecular system:

1. Compute CIS wavefunctions via PySCF (or read them from a SHARC child
   QM interface that produces the same files).
2. Write the wavefunctions in SHARC's native format (`mos`, `dets`,
   `AO_overl.mixed`).
3. Run SHARC's `wfoverlap.x` on them → fragment-state overlap matrices.
4. Compose fragment overlaps into ESD overlaps (`O_det`).
5. Contract with the spin-adaptation matrix `U` and the ECI coefficient
   matrix `C` → ECI state overlaps `S_ij`.

Steps 2–5 are implemented here and tested. Step 1 is provided as a
PySCF-based demo so the pipeline runs without external QM licenses.

## Quick start

    pip install numpy pytest pyscf
    pytest -v                                       # 65 tests
    python examples/build_fragment_overlap.py       # fragment overlap demo
    python examples/full_eci_pipeline_demo.py       # full pipeline demo
    python examples/gap_map_demo.py                 # QM.out section gap map
    python examples/pipeline_report.py              # pipeline stage status

## Requirements

- Python 3.10+, NumPy, pytest
- PySCF for the demo (not required for the module tests)
- A clone of `sharc-md/sharc4` at `~/sharc4` with a compiled
  `wfoverlap.x` (see `BUILDING.md`)

## Layout

    src/
      run_wfoverlap.py      # wrapper around the SHARC wfoverlap binary
      qmout_parser.py       # parser for SHARC's QM.out format + gap map
      esd_to_dets.py        # ESD -> detstring mapper
      fragment_overlap.py   # fragment file writers + wfoverlap invocation
      esd_overlap.py        # ESD overlap product formula
      state_overlap.py      # ECSF contraction + state overlaps
      pipeline_status.py    # pipeline stage reporter
    tests/                  # 65 tests, all passing
    examples/               # runnable demos
    BUILDING.md             # how to build wfoverlap with gfortran
    NOTES.md                # technical writeup of the pipeline and gaps

## Status

| Stage | Status |
|---|---|
| Fragment overlap primitive | ✅ implemented and tested |
| ESD overlap product formula | ✅ implemented and tested |
| ECSF contraction / state overlaps | ✅ implemented and tested |
| Full ECI integration into `SHARC_ECI.py` | ❌ not yet done |
| Energy gradients | ❌ not done |
| Spin-orbit couplings | ❌ not done |

The remaining stages map directly onto the ECI2GAME postdoc deliverables.
