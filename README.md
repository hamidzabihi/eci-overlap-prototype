# eci-overlap-prototype

[![tests](https://github.com/hamidzabihi/eci-overlap-prototype/actions/workflows/test.yml/badge.svg)](https://github.com/hamidzabihi/eci-overlap-prototype/actions/workflows/test.yml)

*Reproduces and extends the ECI method of Piteša et al., [JCTC 2024](https://doi.org/10.1021/acs.jctc.4c00157) (see also [JPCL 2025](https://doi.org/10.1021/acs.jpclett.5c00065)).*

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
    NEXT_MILESTONE.md          # proposed diabatization validation milestone

## Status

| Stage | Status |
|---|---|
| Fragment overlap primitive | ✅ tested |
| ESD overlap product formula | ✅ tested |
| ECSF contraction / state overlaps | ✅ tested |
| Wavefunction phase correction | ✅ tested |
| Diabatization / state tracking | ✅ tested |
| ECI Hamiltonian from fragments | ✅ tested |
| FEM-level ECI (singlet block) | ✅ matches direct CIS (3-frag: 17.6 meV, 4-frag: 19.9 meV MAD) |
| Triplet ECI extension | ⚠️ implemented; diabatic-vs-adiabatic caveat (see NOTES.md) |
| ECIS-level ECI (with GS-LE) | ✅ validated (2-fragment; N-fragment pending) |
| Full SHARC_ECI.py integration | ❌ future work |
| Energy gradients | ⚠️ FD prototype; force-balance validated (see docs/gradients_status.md) |
| Spin-orbit couplings | ❌ future work |

The remaining stages map directly onto the ECI2GAME postdoc deliverables.

## Recent milestone — triplet extension (tag `triplet-extension-v1`)

The ECI Hamiltonian now includes the triplet block alongside the
singlet block, with a 6×6 basis:

    [ GS, S_A, S_B | T_A, T_B, T_A-T_B ]

Key results at 4 Å (ethylene dimer, 6-31G):

| Quantity | Value |
|---|---|
| Triplet Frenkel coupling `V_TT` | −464.80 meV |
| α-only exchange `K` | 17.98 meV |
| T1 (ECI) vs T1 (direct CIS) | 3.1236 eV vs 3.5605 eV (dev −437 meV) |
| Singlet MAD (S1, S2) | 122 meV |
| `2\|V_TT\| = T2 − T1` consistency | holds to < 1 meV |

Separation scan (4–10 Å, singlet Frenkel coupling — the validated
quantity):

| sep (Å) | V_SS (meV) |
|---|---|
| 4.0  | −299.71 |
| 5.0  | −158.85 |
| 6.0  |  −96.25 |
| 8.0  |  −42.62 |
| 10.0 |  −22.31 |

Log-log fit over 4–10 Å: `|V_SS| ~ R^(−2.84)`, approaching the
`R^(−3)` asymptote expected for dipole-dipole coupling.

**Important distinction:** the ECI `V_TT` is a **diabatic** Frenkel
coupling between localized fragment triplet excitations. The direct
CIS T1–T2 splitting is an **adiabatic** splitting of delocalized
combined-system states. These are different observables and agree
only after (a) diabatization of the direct CIS states, and (b)
verification that the fragment triplet and the combined-system
triplet share the same orbital character.

At the RHF/6-31G level, the fragment T1 of ethylene has mixed
`pi -> pi*` and Rydberg character (`<X_singlet, X_triplet> = 0.475`),
so (b) does not hold and the ECI `V_TT` should not be expected to
match the direct CIS T1–T2 splitting. See `NOTES.md` for the full
caveat and `TODO.md` for the follow-up plan (Rydberg projection or
diabatization).

The ECI assumes strongly orthogonal fragments (no inter-fragment
density overlap). For ethylene this holds above ~3.5 Å; the demo
warns below that threshold. See `NOTES.md` for the full validation
summary.

## Known limitations

See `NOTES.md` for a detailed writeup, including:

- The `wfoverlap.x` build issue (gfortran + `-i4` LAPACK causes
  silent numerical corruption for the native-ascii input path).
- The GS-LE coupling is now implemented correctly (V_ne + J − K,
  eq 9 of the ECI paper) and validated on the symmetric ethylene
  dimer: |V_0A| = |V_0B| to machine precision, 12.85 meV at 4 Å.
  The self-consistent EHF embedding (Brillouin's theorem check) is
  not yet implemented. See `docs/gs_le_validation.md`.
- The ECI is only valid for non-overlapping fragments
  (separation ≳ 3.5 Å for the ethylene demo).

## Recent milestone — GS-LE coupling validated

The GS-LE coupling path (previously documented as unreliable) now
implements eq 9 of the ECI paper directly: V_ne + J − K. For a
symmetric ethylene dimer at 4 Å and 6-31G, the two coupling elements
have equal magnitude and opposite sign to machine precision
(12.85 meV each), and the ECI singlet MAD against direct CIS is
~120 meV — consistent with the paper's ECIS + FDA behavior for
weakly coupled systems. Two tests in `tests/test_gs_le_symmetry.py`
guard the fix. Full details in `docs/gs_le_validation.md`.

Open items: EHF embedding (Brillouin's theorem check), strong-coupling
regime (paper's MgG4^2+ case), and multi-state-per-fragment extension.

## Connection to ECI2GAME

The ECI2GAME postdoc project at RBI targets extending the ECI method
with the electronic-structure quantities required for SHARC surface
hopping: energy gradients, wavefunction overlaps, and spin–orbit
couplings. This prototype implements the overlap half of that
pipeline end-to-end:

- **Wavefunction overlaps:** `esd_overlap.py`, `state_overlap.py`, and
  `run_wfoverlap.py` compute the full `S_ij` matrix that SHARC's NAMD
  requires at each step.
- **Diabatization / state tracking:** `diabatization.py` and
  `phase_tracking.py` handle the phase-consistent state assignment
  that overlap-based dynamics needs across a geometry path.
- **ECI Hamiltonian:** `excitonic_hamiltonian.py` reproduces the
  fragment-based Hamiltonian for a two-chromophore test case at the
  FEM level.

Missing pieces (also the ECI2GAME deliverables): gradients, SOCs,
and the full `SHARC_ECI.py` integration. See `NEXT_MILESTONE.md`.

## Related work

- SHARC 4 source: https://github.com/sharc-md/sharc4
- `cis_nto` (CP2K reader PR): https://github.com/marin-sapunar/cis_nto/pull/4
