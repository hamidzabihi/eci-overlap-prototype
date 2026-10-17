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
