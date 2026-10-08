# ECI → SHARC-NAMD interface: status notes

This repository is a **working prototype** of the overlap computation
needed for surface-hopping NAMD with ECI wavefunctions. It implements
the chain from fragment wavefunctions to the final state-overlap matrix
`S_ij` that SHARC consumes.

## What already exists in SHARC

- `lib/ECI.py` — the ECI method. Builds the effective Hamiltonian from
  fragment densities and diagonalizes it. Author: Tomislav Piteša.
- `bin/SHARC_ECI.py` — the SHARC driver for ECI.
- `wfoverlap/source/` — a Fortran engine for Slater-determinant overlaps
  and Dyson orbitals.
- `lib/qmout.py` — the SHARC QM.out format, with ~23 sections.

## What this repository implements

1. **`src/run_wfoverlap.py`** — Python wrapper for `wfoverlap.x`, with a
   regression test against SHARC's own `water_molcas` reference.

2. **`src/qmout_parser.py`** — parser for SHARC's `QM.out` format and a
   gap map showing which sections the ECI interface currently writes.

3. **`src/esd_to_dets.py`** — mapper from ECI Excitonic Slater
   Determinants to the `detstring` format `wfoverlap` reads.

4. **`src/fragment_overlap.py`** — file writers (MO coefficients,
   determinants, AO overlap) and orchestration of the `wfoverlap` call.
   Produces fragment-state overlap matrices.

5. **`src/esd_overlap.py`** — ESD overlap product formula:
   `<ESD_a^A | ESD_b^B> = prod_f <s_f^{a,A} | s_f^{b,B}>`.

6. **`src/state_overlap.py`** — ECSF contraction and state-overlap matrix:
   `S_ij = C_A^H · U_A^T · O_det · U_B · C_B`.

7. **`src/pipeline_status.py`** — reporter showing pipeline stage status.

## What is done

Stages 1–7 of the pipeline are implemented and tested:

- ECI ESD definition (SHARC, pre-existing)
- ESD → detstring mapping
- dets file writer
- wfoverlap engine (SHARC, pre-existing)
- Python wrapper for wfoverlap
- ECSF contraction (`O_det → S_ij`)
- State overlap matrix `S_ij`

## What is missing

Stages 8–11 are not implemented:

- **Stage 8**: Wavefunction phases — sign/phase tracking for the state
  overlap matrix
- **Stage 9**: ECI energy gradients — derivatives of the effective
  Hamiltonian with respect to nuclear coordinates
- **Stage 10**: Non-adiabatic couplings — can be computed from overlaps,
  but the SHARC-side wiring is not done
- **Stage 11**: Spin-orbit couplings — fragment SOC integrals + ECSF
  transformation

These map onto the ECI2GAME postdoc deliverables:
- Wavefunction overlaps → stage 6–8 (mostly done here)
- Energy gradients → stage 9
- Spin-orbit couplings → stage 11

## Format constraints discovered

Two constraints of the SHARC `wfoverlap` interface that any ECI-side
detstring generator must respect:

1. **α/β count preservation.** The Fortran reader (`iomod.f90:
   read_civec_native`) allocates the SSD arrays with the reference
   determinant's α/β counts and does not bounds-check. Every detstring
   in a `dets` file must preserve the α and β counts of the first
   determinant. This is the convention `SHARC_GAUSSIAN.get_dets_from_chk`
   follows, and it must be followed by any ECI detstring generation.

2. **The three-block output.** `wfoverlap` prints three overlap matrices
   in sequence: raw, renormalized, and Löwdin-orthonormalized. A naive
   parser matching `"Overlap matrix"` will also match the
   `"Renormalized overlap matrix"` and `"Orthonormalized overlap matrix"`
   headers. `parse_overlap_matrix` in `src/run_wfoverlap.py` stops at
   the end of the first block.

## Building

See `BUILDING.md` for the gfortran build instructions for `wfoverlap.x`.

## Running

    pytest -v                                      # 65 tests
    python examples/build_fragment_overlap.py      # fragment overlap demo
    python examples/full_eci_pipeline_demo.py      # full pipeline demo
    python examples/gap_map_demo.py                # QM.out gap map
    python examples/pipeline_report.py             # pipeline stage status
