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
