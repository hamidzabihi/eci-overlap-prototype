# eci-overlap-prototype

Prototype bridge between the ECI module of SHARC and the existing
`wfoverlap` Fortran engine, aimed at enabling surface-hopping
nonadiabatic dynamics with ECI wavefunctions.

## Status

Phase 1 (wfoverlap wrapper + regression test) is complete.

## Tests

    pytest -v

Requires a clone of `sharc-md/sharc4` at `~/sharc4` with a compiled
`wfoverlap.x` (see BUILDING.md).
