#!/usr/bin/env bash
set -euo pipefail
: "${ECI_PROTO:?Set ECI_PROTO}"

GAME="$ECI_PROTO/bin/SHARC_ECI_GAME.py"
[[ -f "$GAME" ]] || { echo "missing $GAME" >&2; exit 1; }
python -c "import ast; ast.parse(open('$GAME').read())"
echo "ok: $GAME parses"

mkdir -p "$ECI_PROTO/docs"

cat > "$ECI_PROTO/docs/ECI2GAME_status.md" <<'MDEOF'
# ECI2GAME - Phase 1 status

| Deliverable | Status | Evidence | Next step |
|---|---|---|---|
| Wavefunction overlaps (S_ij) | prototype + integrated | src/esd_overlap.py, src/state_overlap.py, bin/SHARC_ECI_GAME.py | - |
| Phase tracking | done | src/phase_tracking.py (matches utils.phase_correction_cmplx) | - |
| Diabatization / state tracking | done | src/diabatization.py | - |
| ECI Hamiltonian (FEM) | done (13-20 meV MAD, 2-4 frag) | README status table | - |
| ECIS with GS-LE | done (2-frag validated) | docs/gs_le_validation.md | N-fragment pending |
| Energy gradients | FD prototype only | docs/gradients_status.md | Analytic - multi-week |
| Spin-orbit couplings | not started | - | Multi-week |
| Full SHARC_ECI integration | Phase-1 subclass | bin/SHARC_ECI_GAME.py | Phase-2 upstream PR |
| Triplet ECI | experimental | README, NOTES.md | Rydberg projection |

**Phase-1 scope boundary:** overlap is only defined when the ECI basis is
identical between consecutive steps. Basis-changing overlap raises
`NotImplementedError` in `ECI.calculate_overlap_game`.

**Environment:** `SHARC_ECI_GAME.py` requires `$SHARC/bin` and `$SHARC/lib`
on `PYTHONPATH` because it reuses `SHARC_INTERFACE`, `SHARC_HYBRID`,
`SHARC_ECI`, and `ECI` from the sharc4 installation. Zero sharc4 files are
modified; the extension is a subclass of `SHARC_ECI` plus monkeypatches of
`ECI.ECI.__init__` and `ECI.ECI.run`.
MDEOF

cat > "$ECI_PROTO/docs/upstream_pr.md" <<'MDEOF'
# Phase-2 upstream PR to sharc-md/sharc4

Minimal diff once the prototype is validated.

## lib/ECI.py

- `ECI.ECI.__init__(..., previous=None)`: store `self.previous`, add
  `self.overlap = {m: None for m in job.multiplicities}`.
- `ECI.ECI.run()`: call `self.calculate_overlap()` between
  `calculate_eigenstates()` and `calculate_properties()` if
  `self.previous is not None`.
- New methods `ECI.ECI.calculate_overlap(self)` and
  `ECI.ECI.get_state(self)` - see `bin/SHARC_ECI_GAME.py`.

## bin/SHARC_ECI.py

- `all_features` -> `{"h", "dm", "mol", "point_charges", "overlap", "phases"}`
- `read_requests`: propagate `overlap` to SSC children, imply
  `overlap` from `phases`.
- `run()`: load/save `eci_wf.<Z>.<step>.pkl` from `savedir`.
- `getQMout()`: populate `QMout["overlap"]` from `ECIjobs[Z].overlap`
  and `QMout["phases"]` via `utils.phase_correction_cmplx(..., flag="simple")`.

Plus `tests/interface/test_eci_overlap.py` following the existing
`test_interface.py` structure.

Total: ~30 lines + test.
MDEOF

echo "wrote docs/ECI2GAME_status.md and docs/upstream_pr.md"
