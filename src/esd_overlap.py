"""
Compute ESD overlaps from fragment-state overlaps.

An Excitonic Slater Determinant (ESD) in ECI is a product of per-fragment
electronic states:

    |ESD> = |s_A> (x) |s_B> (x) ...

The overlap between two ESDs at different geometries is the product of
per-fragment overlaps:

    <ESD_a^A | ESD_b^B> = prod_f <s_f^{a,A} | s_f^{b,B}>

This module implements that product formula, taking a dictionary of
fragment-level state overlaps as input and returning the ESD-ESD
overlap matrix.

Design notes
------------
The fragment-state overlaps O^(f)[s_A, s_B] must be computed externally
(e.g. via src.fragment_overlap for a fragment at two geometries). This
module is agnostic to how they were obtained.

An ESD is represented here as a list of (fragment_label, state_index)
tuples, one per fragment. The state_index is the ordinal within the
fragment's state list for its charge state.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Sequence
import numpy as np


# ---------------------------------------------------------------- ESD model

@dataclass(frozen=True)
class ESDSpec:
    """Specification of an ESD as a list of (fragment, state_index) pairs.

    Attributes
    ----------
    site_states : tuple of (fragment_label, state_index)
        Each entry identifies one fragment and the index of the electronic
        state it is in for this ESD. The order of the tuple defines the
        canonical ordering of fragments in the ESD.
    """
    site_states: tuple

    def __post_init__(self):
        if not self.site_states:
            raise ValueError("ESD must have at least one fragment")

    @property
    def fragments(self) -> tuple:
        return tuple(f for f, _ in self.site_states)

    def state_of(self, fragment: str) -> int:
        for f, s in self.site_states:
            if f == fragment:
                return s
        raise KeyError(f"fragment {fragment!r} not in ESD")

    def __repr__(self) -> str:
        parts = [f"{f}[{s}]" for f, s in self.site_states]
        return "ESD(" + ", ".join(parts) + ")"


# ---------------------------------------------------------------- the product

def esd_overlap(
    esd_a: ESDSpec,
    esd_b: ESDSpec,
    fragment_overlaps: dict,
) -> complex:
    """Compute <ESD_a^A | ESD_b^B> as the product of fragment-state overlaps.

    Parameters
    ----------
    esd_a, esd_b : ESDSpec
        The two ESDs.
    fragment_overlaps : dict
        Maps (fragment_label, state_a_index, state_b_index) -> overlap value.
        For example: fragment_overlaps[("BD", 0, 0)] = 0.997.
        Must contain an entry for every fragment in esd_a (and esd_b).

    Returns
    -------
    complex: the ESD overlap.
    """
    if esd_a.fragments != esd_b.fragments:
        raise ValueError(
            f"ESD fragment sets differ: {esd_a.fragments} vs {esd_b.fragments}"
        )

    result = complex(1.0, 0.0)
    for frag, s_a in esd_a.site_states:
        s_b = esd_b.state_of(frag)
        key = (frag, s_a, s_b)
        if key not in fragment_overlaps:
            raise KeyError(f"missing fragment overlap {key}")
        result *= fragment_overlaps[key]
    return result


# ---------------------------------------------------------------- the matrix

def esd_overlap_matrix(
    esds_a: Sequence[ESDSpec],
    esds_b: Sequence[ESDSpec],
    fragment_overlaps: dict,
) -> np.ndarray:
    """Compute the full (nA x nB) matrix of ESD overlaps.

    Parameters
    ----------
    esds_a, esds_b : sequence of ESDSpec
    fragment_overlaps : dict
        See esd_overlap.

    Returns
    -------
    (nA, nB) complex matrix O_det where O_det[i, j] = <ESD_i^A | ESD_j^B>.
    """
    nA = len(esds_a)
    nB = len(esds_b)
    O = np.zeros((nA, nB), dtype=complex)
    for i, a in enumerate(esds_a):
        for j, b in enumerate(esds_b):
            O[i, j] = esd_overlap(a, b, fragment_overlaps)
    return O


# ---------------------------------------------------------------- helpers

def fragment_overlaps_from_matrix(
    fragment_label: str,
    overlap_matrix: np.ndarray,
) -> dict:
    """Convert a fragment overlap matrix (nA, nB) into the dict form used here.

    Parameters
    ----------
    fragment_label : str
        The fragment label (e.g. "BD").
    overlap_matrix : (nA, nB) array
        Matrix of state overlaps for this fragment.

    Returns
    -------
    dict mapping (fragment_label, i, j) -> overlap_matrix[i, j]
    """
    nA, nB = overlap_matrix.shape
    return {
        (fragment_label, i, j): complex(overlap_matrix[i, j])
        for i in range(nA)
        for j in range(nB)
    }


def merge_fragment_overlaps(*dicts: dict) -> dict:
    """Merge several fragment-overlap dictionaries into one.

    Useful when each fragment's overlaps are computed separately.
    """
    merged = {}
    for d in dicts:
        merged.update(d)
    return merged
