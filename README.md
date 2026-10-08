# eci-overlap-prototype

A prototype exploring the path from **ECI** (Excitonic Configuration
Interaction) to **SHARC's surface-hopping NAMD**, by building the
scaffolding around the interface gaps.

## What this is

Not a full implementation of the ECI-NAMD interface. It is a
**status analysis** backed by working code: it traces the pipeline
from an ECI wavefunction to the SHARC `QM.out` format, verifies each
existing piece against real SHARC data, and enumerates exactly what
is missing.

The missing pieces map 1:1 onto the **ECI2GAME** postdoc deliverables
(energy gradients, wavefunction overlaps, spin-orbit couplings).

## Quick start

    pytest -v                              # 29 tests
    python examples/gap_map_demo.py        # QM.out section gap map
    python examples/pipeline_report.py     # pipeline stage status

## Requirements

- Python 3.10+, NumPy, pytest
- A clone of `sharc-md/sharc4` at `~/sharc4` with a compiled
  `wfoverlap.x` (see `BUILDING.md`)

## Layout

    src/
      run_wfoverlap.py      # Python wrapper around the SHARC wfoverlap binary
      qmout_parser.py       # parser for SHARC's QM.out + gap map
      esd_to_dets.py        # ECI ESD -> wfoverlap detstring mapper
      pipeline_status.py    # pipeline stage reporter
    tests/                  # 29 tests, all passing
    examples/               # runnable demos
    BUILDING.md             # how to build wfoverlap with gfortran
    NOTES.md                # technical writeup
