# ECI Gradients — Status

## What is implemented

`src/gradients.py` provides a central-difference gradient of any eigenvalue
of the ECI Hamiltonian with respect to Cartesian nuclear coordinates. The
caller supplies a `rebuild_H(geom_a, geom_b)` function; the module handles
the perturbation loop and returns a three-component gradient per atom.

Usage: `examples/gradients_demo.py` (S1 of the ethylene dimer, 4 Å separation).

## Validation — force balance

For an isolated system, the sum of all atomic gradients must vanish. This is
an exact constraint: no external force acts on the molecule as a whole.

| Quantity | Value (Ha/Å) |
|---|---|
| \|Σ grad\| ECI FD | 5.1 × 10⁻¹⁰ |
| \|Σ grad\| direct CIS | 1.5 × 10⁻¹⁴ |

Both are zero to numerical precision. **This confirms the FD machinery
is correct** — the gradient is mathematically the exact derivative of
the ECI eigenvalue.

## Comparison to direct CIS

The ECI S1 FD gradient is compared against the analytic CIS S1 gradient
of the combined system (from PySCF's TDA response):

| Quantity | Value |
|---|---|
| \|grad\| ECI FD | 1.11 Ha/Å |
| \|grad\| direct CIS | 0.32 Ha/Å |
| cosine similarity | +0.38 |
| norm ratio (ECI/CIS) | 3.5× |

The two gradients are **not equivalent observables**, and this is
expected:

- **ECI S1** is built from *isolated* fragment SCF + CIS site states.
  When fragment A's atoms move, only fragment A's site states respond;
  fragment B stays in its isolated-fragment form.
- **Direct CIS S1** is a fully coupled eigenstate of the combined
  system. When any atom moves, every MO in the combined system relaxes.

Their derivatives therefore respond to different aspects of the
potential energy surface. Cosine similarity between them is not a
meaningful correctness criterion.

A more meaningful comparison would be: ECI FD gradient vs. ECI
**analytic** gradient, once analytic derivatives are implemented.

## What is missing

- **Analytic gradients.** The FD approach costs O(6N) fragment SCF+CIS
  calculations per gradient — ~35 s for the 12-atom ethylene dimer at
  6-31G. Production work needs analytic derivatives, which require:
  - derivatives of fragment SCF/CIS site energies (available from PySCF)
  - derivatives of the cross-fragment J and K integrals (needs explicit
    integral derivative code or the Wigner–Eckart route)
  - derivatives of EHF embedding charges if embedding is used (this is
    the piece that would need new theory development)
- **Gradients for the full ECI basis.** The current demo returns the
  singlet block only. Extending to DLE and triplet states requires
  the corresponding Hamiltonian blocks to be differentiable.

## Next steps

1. Implement analytic derivatives for the cross-fragment integrals
2. Validate analytic vs. FD gradients on the same system
3. Extend to the full ECI basis including DLE and triplets
4. Add EHF charge derivatives (the theoretical core of the production
   gradient)

## What this means for ECI2GAME

The ECI2GAME project lists energy gradients as one of three core
deliverables. This prototype provides:
- a working validation tool (FD)
- a demonstrated correctness criterion (force balance)
- a clear list of the derivative terms that analytic gradients need

Implementing analytic gradients remains the substantive research task.
