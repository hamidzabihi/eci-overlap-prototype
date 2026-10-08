# ECI -> SHARC-NAMD interface: status notes

This repository is a **prototype** exploring the path from the ECI
(Excitonic Configuration Interaction) method to SHARC's non-adiabatic
molecular dynamics (NAMD) machinery. It does **not** implement the full
ECI-NAMD interface — it lays out the scaffolding around the gaps, so
that the gaps themselves become unambiguous.

## What already exists in SHARC

- `lib/ECI.py` — the ECI method. Builds the effective Hamiltonian from
  fragment densities and diagonalizes it. Author: Tomislav Piteša.
- `bin/SHARC_ECI.py` — the SHARC driver for ECI. Orchestrates EHF and
  site-state child jobs, then calls `ECI.ECI`.
- `wfoverlap/source/` — a Fortran engine that computes Slater-determinant
  overlaps and Dyson orbitals. Mature, tested, independent of ECI.
- `lib/qmout.py` — the SHARC `QM.out` format, defining ~23 sections
  (`! 0` through `! 999`).

## What this repository implements

Three small modules, each tested against real SHARC data or the real
SHARC format:

1. **`src/run_wfoverlap.py`** — a Python wrapper for the `wfoverlap.x`
   binary. Four regression tests validate it against SHARC's own
   `water_molcas` reference output.

2. **`src/qmout_parser.py`** — a parser for the SHARC `QM.out` format.
   Six tests validate it against the real BODIPY-dimer ECI output
   (`examples/SHARC_ECI/QM.out`). The parser's `gap_report()` shows
   which sections the ECI interface currently writes.

3. **`src/esd_to_dets.py`** — a mapper from ECI Excitonic Slater
   Determinants to the `detstring` format `wfoverlap` reads. Thirteen
   tests validate the mapping, electron bookkeeping, and file round-trip.

4. **`src/pipeline_status.py`** — a reporter that walks through the
   ECI -> wfoverlap -> SHARC-NAMD pipeline and marks each stage as
   present, done, or missing.

## What is missing (and why it matters)

Running `python examples/pipeline_report.py` shows that stages 6–11
of the pipeline are not implemented anywhere:

- Stage 6: ECSF contraction (`O_det -> S_ij`)
- Stage 7: State overlap matrix `S_ij`
- Stage 8: Wavefunction phases
- Stage 9: ECI energy gradients
- Stage 10: Non-adiabatic couplings
- Stage 11: Spin-orbit couplings

Running `python examples/gap_map_demo.py` confirms the same picture
from the SHARC side: the ECI interface writes only `QM.out` sections
0, 1, 2, and 8. Sections 3, 5, 6, 7, and 13 are absent — and those
are precisely the sections SHARC reads for surface-hopping NAMD.

## Mapping to the ECI2GAME postdoc deliverables

The ECI2GAME posting asks for work on:

1. **ECI energy gradients** — pipeline stage 9, QM.out section 3
2. **Wavefunction overlaps** — pipeline stages 6–8, QM.out sections 6–7
3. **Spin-orbit couplings** — pipeline stage 11, QM.out section 13

The gap between the current code and the postdoc's deliverables is
therefore **explicit, verifiable, and small** — not a research unknown.

## Building

See `BUILDING.md` for the gfortran build instructions for `wfoverlap.x`.

## Running

    pytest -v                              # 29 tests
    python examples/gap_map_demo.py        # QM.out section gap map
    python examples/pipeline_report.py     # pipeline stage status
