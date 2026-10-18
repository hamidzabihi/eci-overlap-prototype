## Real SHARC QM.out validation (Oct 18)

### Corpus — 8/8 files parse
All eight test files in `~/sharc4/tests/interface/` parse through
`src/qmout_parser.py` with correct Hamiltonian shapes:

| File | natom | states | `H` shape | Special sections |
|---|---|---|---|---|
| QM.out | 6 | [4] | (4,4) | !3 grads, !5 NACs, !30/!31 |
| qmout_dens.out | 3 | [2] | (2,2) | !22, !24 |
| qmout_ion.out | 16 | [6,10,8] | (50,50) | !20 property matrices |
| qmout_nacdr.out | 3 | [4] | (4,4) | !5 only |
| qmout_notes.out | 4 | [4] | (4,4) | !999 |
| qmout_overlap.out | 4 | [4] | (4,4) | !6 overlap, !7 phases |
| qmout_theodore.out | 6 | [2] | (2,2) | !3, !21 |
| QM1.out (MNDO) | 6 | [4] | (4,4) | !3, !5, !30, !31 |

Coverage: natom 3–16, nmstates 2–50, singlets/doublets/triplets,
with/without point charges. All 8 files also print cleanly through
SHARC's own `QMout_print.py`.

### SHARC reader note
`SHARC_QMOUT.py` exposes `SHARC_QMOUT` as a *writer* class only — no
module-level read function. The canonical reader is `qmout.QMout`
(used by `QMout_print.py`). In SHARC 4, `qm.prop0d/prop1d/prop2d`
stay empty after construction; a small raw parser for `! 6` / `! 7`
was written instead (see `scripts/read_overlap.py`).

### Ground-truth integration with the prototype — CLOSED
`qmout_overlap.out` parsed:

- `! 6` overlap  S:  (4, 4) complex128
- `! 7` phases   :  (4,)   complex128, all values `1+0j`

Fed into prototype modules:

| Call | Result |
|---|---|
| `so.orthonormality_check(S)` | `ok=False`, residual = `3.271e-03` — correct: S is a real overlap between slightly displaced geometries, not exactly I |
| `pt.simple_phase_correction(S)` | phases = `[1,1,1,1]`; diag real-positive after correction |
| `pt.apply_phases(S, phases)` | `|dS|_max = 0.0` |
| `pt.max_diagonal_imag(S)` | `0.0` |
| `pt.min_diagonal_magnitude(S)` | `1.000000` |

Unitarity check on real S:
- `|S - I|_max = 3.271e-03`
- `|S*S^H - I|_max = 6.539e-03`

These residuals are the physical geometry drift in the test file,
not errors. The prototype reads, checks, and phase-corrects the real
SHARC overlap data correctly.

### Local evidence summary
- 9/9 capability validator
- 93 pytest tests
- 8/8 SHARC `QM.out` test files parsed by `src/qmout_parser.py`
- 8/8 validated through SHARC's `QMout_print.py`
- Ground-truth `!6`/`!7` from `qmout_overlap.out` fed into
  `src/state_overlap.py` and `src/phase_tracking.py`
