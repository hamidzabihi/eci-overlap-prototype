"""
Map ECI Excitonic Slater Determinants (ESDs) to wfoverlap detstrings.

Background
----------
An ESD in ECI is a product of per-fragment electronic states:

    |ESD> = |s_A> (x) |s_B> (x) ...

Each fragment state |s_X> is characterized by:
    - a charge Z_X
    - a spin multiplicity 2*S_X (0 = singlet, 1 = doublet, 2 = triplet, ...)
    - an ordinal N_X within its (Z, S) manifold
    - a set of occupied MOs in the fragment's MO basis

The full-system detstring that wfoverlap reads has length nMO_total,
where nMO_total = sum over fragments of nMO_fragment. Each fragment's
MOs form a contiguous block. Within a block, each MO is labeled:

    'd' = doubly occupied
    'a' = alpha electron
    'b' = beta electron
    'e' = empty

This module provides the mapping given explicit occupation information.
It does NOT extract occupations from a fragment's CI vector — that
would require natural orbital analysis and is a research-grade task
for the ECI2GAME postdoc itself.

dets file format (from SHARC's read_civfl.py, write_det_file):
    nstate nMO nSD
    detstring  coeff1  coeff2  ...  coeff_nstate
    ...
Coefficients are formatted as ' % 14.10f' (i.e., 17 chars wide).
"""

from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence
import numpy as np


@dataclass
class FragmentOccupation:
    """Occupation pattern of one fragment in one electronic state.

    Attributes
    ----------
    label     : fragment label (e.g., "BD", "BA")
    Z         : total charge of the fragment state
    S         : two times the spin quantum number (0 for singlet, 1 for doublet, ...)
    N         : ordinal within the (Z, S) manifold (1-based)
    n_mo      : total number of MOs in this fragment's basis
    alpha_occ : 0-based indices of occupied alpha MOs
    beta_occ  : 0-based indices of occupied beta MOs
    """
    label: str
    Z: int
    S: int
    N: int
    n_mo: int
    alpha_occ: list[int] = field(default_factory=list)
    beta_occ: list[int] = field(default_factory=list)

    def __post_init__(self):
        for i in self.alpha_occ:
            if not 0 <= i < self.n_mo:
                raise ValueError(
                    f"alpha_occ index {i} out of range for {self.label} (n_mo={self.n_mo})"
                )
        for i in self.beta_occ:
            if not 0 <= i < self.n_mo:
                raise ValueError(
                    f"beta_occ index {i} out of range for {self.label} (n_mo={self.n_mo})"
                )

    @property
    def n_alpha(self) -> int:
        return len(self.alpha_occ)

    @property
    def n_beta(self) -> int:
        return len(self.beta_occ)


def fragment_to_occupation(frag: FragmentOccupation) -> list[str]:
    """Convert one fragment's occupation into per-MO symbols of length n_mo.

    Rules
    -----
    - an MO in alpha_occ but not beta_occ -> 'a'
    - an MO in beta_occ but not alpha_occ -> 'b'
    - an MO in both                      -> 'd'
    - an MO in neither                   -> 'e'

    Returns
    -------
    list of single-character strings, one per MO in this fragment's block.
    """
    occ = ["e"] * frag.n_mo
    for i in frag.alpha_occ:
        occ[i] = "a"
    for i in frag.beta_occ:
        if occ[i] == "a":
            occ[i] = "d"
        elif occ[i] == "e":
            occ[i] = "b"
        else:
            # should be unreachable given __post_init__ validation
            raise ValueError(
                f"inconsistent occupation in {frag.label} at MO {i}"
            )
    return occ


def esd_to_detstring(fragments: Sequence[FragmentOccupation]) -> str:
    """Concatenate per-fragment occupations into a single detstring.

    The fragment ordering defines the MO ordering in the detstring.
    All fragments must be listed in the same order as the ESD's
    site_states dict iteration.
    """
    return "".join(s for frag in fragments for s in fragment_to_occupation(frag))


def total_electrons(fragments: Sequence[FragmentOccupation]) -> tuple[int, int]:
    """Return (n_alpha, n_beta) for the concatenated detstring."""
    n_alpha = sum(f.n_alpha for f in fragments)
    n_beta = sum(f.n_beta for f in fragments)
    return n_alpha, n_beta


def total_charge(fragments: Sequence[FragmentOccupation]) -> int:
    """Sum of fragment charges."""
    return sum(f.Z for f in fragments)


def total_nmo(fragments: Sequence[FragmentOccupation]) -> int:
    """Total number of MOs across all fragments."""
    return sum(f.n_mo for f in fragments)


def total_spin_S(fragments: Sequence[FragmentOccupation]) -> int:
    """Sum of 2*S over fragments. This is an upper bound on the total
    spin; the actual total S depends on how fragment spins couple."""
    return sum(f.S for f in fragments)


# ---------------------------------------------------------------- dets writer

DETS_COEFF_FMT = " % 14.10f"


def write_dets_file(
    path: Path | str,
    detstrings: Sequence[str],
    coefficients: np.ndarray,
    nmo: int,
) -> None:
    """Write a dets file in the format wfoverlap expects.

    Parameters
    ----------
    path         : output file path
    detstrings   : list of determinant strings (each of length nmo)
    coefficients : array of shape (nstate, nSD) — coefficients[i, d] is the
                   coefficient of determinant d in state i. Note the transpose
                   relative to some conventions.
    nmo          : number of MOs (must match len(detstring))

    The written file has the form:
        nstate nmo nSD
        detstring  c0  c1  ...  c_{nstate-1}
        ...
    """
    path = Path(path)
    if coefficients.ndim != 2:
        raise ValueError("coefficients must be a 2D array (nstate, nSD)")
    nstate, nsd = coefficients.shape
    if nsd != len(detstrings):
        raise ValueError(
            f"coefficients has {nsd} determinants, detstrings has {len(detstrings)}"
        )
    for i, d in enumerate(detstrings):
        if len(d) != nmo:
            raise ValueError(
                f"detstring {i} has length {len(d)}, expected {nmo}"
            )

    lines = [f"{nstate} {nmo} {nsd}"]
    for d, det in enumerate(detstrings):
        coeffs = "".join(DETS_COEFF_FMT % coefficients[s, d] for s in range(nstate))
        lines.append(det + coeffs)
    path.write_text("\n".join(lines) + "\n")


def read_dets_file(path: Path | str) -> tuple[int, int, list[str], np.ndarray]:
    """Read a dets file. Returns (nstate, nmo, detstrings, coefficients)."""
    path = Path(path)
    lines = path.read_text().splitlines()
    header = lines[0].split()
    nstate, nmo, nsd = int(header[0]), int(header[1]), int(header[2])
    detstrings: list[str] = []
    coeffs = np.zeros((nstate, nsd))
    for d in range(nsd):
        line = lines[1 + d]
        detstrings.append(line[:nmo])
        rest = line[nmo:]
        # Parse coefficients by whitespace, matching Fortran's list-directed
        # read (READ(coefstring,*)). This is robust to any field width.
        tokens = rest.split()
        if len(tokens) != nstate:
            raise ValueError(
                f"determinant {d}: expected {nstate} coefficients, "
                f"found {len(tokens)}"
            )
        for s in range(nstate):
            coeffs[s, d] = float(tokens[s])
    return nstate, nmo, detstrings, coeffs


def detstring_to_symbol_counts(detstring: str) -> dict[str, int]:
    """Count the number of each symbol in a detstring. Useful for tests."""
    return {
        "d": detstring.count("d"),
        "a": detstring.count("a"),
        "b": detstring.count("b"),
        "e": detstring.count("e"),
    }
