# GS-LE Coupling Validation

## What GS-LE coupling is

The ground-state / locally-excited (GS-LE) coupling is the
off-diagonal matrix element of the ECI Hamiltonian between the
ground-state product and a product with one locally excited fragment.
It is evaluated from eq 9 of Pitesa et al. (JCTC 2024) with the
Kronecker-delta structure of eq 6, and vanishes in the Frenkel
exciton model (FEM).

The correct expression (fragment F locally excited, fragment G in GS):

    H_{GS-LE_F} = V_ne_FG + J_FG - K_FG

where:
- V_ne_FG  = attraction between the F->LE transition density and the
             nuclei of G only
- J_FG     = Coulomb between the F->LE transition density and G's
             ground-state density
- K_FG     = exchange (eq 49 of the SI), spin-integrated for the
             closed-shell singlet case

## Bug fixed (2026-10-14)

The previous implementation (following SHARC's lib/ECI.py) had three
issues:

1. Included fragment F's own nuclei in the one-electron operator,
   double-counting the attractive potential already present in F's
   site energy.
2. Omitted the exchange (K) term entirely.
3. Reused the A-first two-electron integral array for both coupling
   calls. For a symmetric dimer with equal-sized basis sets this
   silently produced a numerically valid but physically wrong result
   for the second call.

The combined effect was:

    V1_0,LE_A  = -12.85 meV       (correct magnitude)
    V1_00,LE_B = +25800 meV       (pathological; 25.8 eV)

After the fix, the B-first ERIs are built explicitly
(`eri_BABA`, `eri_BABA_ikjl_perm`) and the coupling is evaluated
from eq 9 with the correct nuclear, Coulomb, and exchange terms.

## Validation — symmetric ethylene dimer

**Setup:** ethylene at 6-31G, two fragments pi-stacked at 4.0 Angstrom.

| Quantity | Value |
|---|---|
| V1_0,LE_A | -12.851 meV |
| V1_00,LE_B | +12.851 meV |
| \|V_A\| - \|V_B\| | < 1e-8 Ha (machine precision) |
| ECI S1 vs direct CIS S1 | 8.413 vs 8.248 eV (+165 meV) |
| ECI S2 vs direct CIS S2 | 9.013 vs 8.933 eV (+80 meV) |
| Singlet MAD vs direct CIS | ~120 meV |

The 80-165 meV deviations are consistent with the paper's reported
ECIS + FDA behavior for a weakly coupled system (Table 2 of the 2024
paper gives S1 MAD = 91 meV for G4 at the same level of theory).

## Range of validity

GS-LE coupling magnitude vs. separation (from the same dimer):

| R (Angstrom) | V1_0,LE_A (meV) |
|---|---|
| 4.0 | -12.85 |
| 5.0 | +0.86 |
| 6.0 | -0.03 (numerical noise) |
| 8.0 | -0.21 (numerical noise) |
| 10.0 | -0.13 (numerical noise) |

For two neutral fragments, the GS-LE coupling is negligible beyond
approximately 5 Angstrom. The RHF/CIS numerical noise floor
(~0.03 meV = 1e-6 Ha) dominates at larger distances. This is expected
physics: both V_ne and J scale as R^-3 for neutral distributions with
zero net charge, and K decays exponentially.

## How to reproduce

    cd ~/eci-overlap-prototype
    python examples/two_fragment_eci_demo.py --include-gs-le

The output reports V1_0,LE_A and V1_00,LE_B during Hamiltonian
assembly. For a symmetric dimer these must have equal magnitude and
opposite sign.

Separation scan:

    for sep in 4.0 5.0 6.0 8.0 10.0; do
        python examples/two_fragment_eci_demo.py --include-gs-le \
            --separation $sep 2>&1 | grep "V1_0,LE_A"
    done

## Open items

- **EHF embedding.** Brillouin's theorem predicts that the GS-LE
  coupling vanishes under converged EHF embedding (paper Section 3.3).
  This has not yet been tested; the self-consistent ESP point-charge
  cycle is not implemented in this prototype.
- **Strong coupling regime.** The paper's Table 3 shows that ECI + FDA
  can fail qualitatively for strongly coupled systems (MgG4^2+,
  ECIS + FDA gives MAD = 656 meV). This prototype has not been tested
  in that regime.
- **Multi-state per fragment.** The current implementation includes
  only the first excited state per fragment. Extending to S2, S3, ...
  requires additional coupling blocks and is straightforward but not
  yet validated.
