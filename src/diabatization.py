"""
State assignment / diabatization by maximum overlap.

Given adiabatic states at a reference geometry R0 and at a displaced
geometry R1, assign each state at R1 to the reference state with the
largest |overlap|. This is the "reference-geometry" tracking algorithm
described in:

    Sapunar, Pitesa, Davidovic, Doslic,
    J. Chem. Theory Comput. 2019, 15, 3461-3469 (section 3.3)

The algorithm is the basis for:
  - state tracking along NAMD trajectories
  - decomposing absorption spectra into diabatic character
  - assigning excited states across different electronic structure
    methods

This module provides:
  - `assign_states(S)` - greedy one-to-one assignment from an overlap matrix
  - `track_states(S_list)` - assignment across a path of geometries
  - `assignment_table(...)` - human-readable summary
"""

from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np


# ---------------------------------------------------------------- single step

@dataclass
class Assignment:
    """One-to-one assignment of new states to reference states.

    Attributes
    ----------
    mapping : dict[int, int]
        mapping[j] = i means new state j is assigned to reference state i
        (both 0-based). -1 means unassigned.
    overlaps : dict[int, float]
        overlaps[j] = |S[i, j]| for the assigned pair (i = mapping[j])
    reference_indices : list[int]
        sorted list of reference indices i that were assigned
    new_indices : list[int]
        sorted list of new indices j that were assigned
    unassigned_new : list[int]
        new indices that could not be assigned (all reference states taken)
    """
    mapping: dict
    overlaps: dict
    reference_indices: list
    new_indices: list
    unassigned_new: list = field(default_factory=list)


def assign_states(S: np.ndarray, threshold: float = 0.1) -> Assignment:
    """Greedy one-to-one assignment of new states (columns) to reference
    states (rows) based on the largest |overlap|.

    Parameters
    ----------
    S : (n_ref, n_new) complex overlap matrix
    threshold : minimum |S_ij| for an assignment to be accepted. New
        states whose best available reference is below this threshold
        are left unassigned (mapped to -1), signalling that they have
        no clear reference character. Default 0.1, which is small enough
        to catch nearly-orthogonal states but large enough to reject
        numerical noise.

    Returns
    -------
    Assignment object

    Notes
    -----
    Greedy: iteratively pick the global maximum |S_ij| among unassigned
    pairs. This is optimal when no two new states prefer the same
    reference state (the usual case for adiabatic tracking). For
    pathological cases, use scipy.optimize.linear_sum_assignment.
    """
    S = np.asarray(S)
    if S.ndim != 2:
        raise ValueError("S must be 2D")
    n_ref, n_new = S.shape

    # Collect all candidate (|S|, i, j) triples above threshold
    triples = []
    for i in range(n_ref):
        for j in range(n_new):
            val = abs(S[i, j])
            if val > threshold:
                triples.append((val, i, j))
    triples.sort(reverse=True)

    mapping = {j: -1 for j in range(n_new)}
    overlaps = {j: 0.0 for j in range(n_new)}
    used_ref = set()

    for val, i, j in triples:
        if mapping[j] != -1:
            continue
        if i in used_ref:
            continue
        mapping[j] = i
        overlaps[j] = val
        used_ref.add(i)

    reference_indices = sorted(set(mapping.values()) - {-1})
    new_indices = sorted([j for j, i in mapping.items() if i != -1])
    unassigned_new = sorted([j for j, i in mapping.items() if i == -1])

    return Assignment(
        mapping=mapping,
        overlaps=overlaps,
        reference_indices=reference_indices,
        new_indices=new_indices,
        unassigned_new=unassigned_new,
    )


# ---------------------------------------------------------------- path

@dataclass
class StateTrack:
    """Track of one reference state across a path of geometries."""
    reference_index: int
    # list of (step_index, matched_new_index, overlap)
    steps: list = field(default_factory=list)


def track_states(S_list: list[np.ndarray]) -> list[Assignment]:
    """Assign states across a path.

    Parameters
    ----------
    S_list : list of (n_ref, n_new) overlap matrices, one per step
        Each matrix compares the SAME reference geometry R0 to the
        geometry at step k.

    Returns
    -------
    list of Assignment, one per step
    """
    return [assign_states(S) for S in S_list]


def track_reference_state(
    assignments: list[Assignment],
    reference_index: int,
) -> StateTrack:
    """Extract the trajectory of a single reference state across steps."""
    track = StateTrack(reference_index=reference_index)
    for step, a in enumerate(assignments):
        for j, i in a.mapping.items():
            if i == reference_index:
                track.steps.append((step, j, a.overlaps[j]))
                break
    return track


# ---------------------------------------------------------------- reporting

def assignment_table(
    assignment: Assignment,
    reference_labels: list[str] | None = None,
    new_labels: list[str] | None = None,
) -> str:
    """Human-readable assignment table.

    Parameters
    ----------
    assignment : Assignment
    reference_labels : optional labels for reference states
    new_labels : optional labels for new states

    Returns
    -------
    multi-line string
    """
    lines = []
    lines.append(f"{'new':<8} {'->':<4} {'ref':<8} {'|overlap|':>12}")
    lines.append("-" * 36)
    for j in sorted(assignment.mapping):
        i = assignment.mapping[j]
        ov = assignment.overlaps[j]
        new_lbl = new_labels[j] if new_labels else f"S{j+1}"
        if i == -1:
            lines.append(f"{new_lbl:<8} {'-':<4} {'?':<8} {'-':>12}")
        else:
            ref_lbl = reference_labels[i] if reference_labels else f"S{i+1}"
            lines.append(f"{new_lbl:<8} {'->':<4} {ref_lbl:<8} {ov:>12.6f}")
    if assignment.unassigned_new:
        lines.append("")
        lines.append(f"Unassigned new states: "
                     + ", ".join(str(j+1) for j in assignment.unassigned_new))
    return "\n".join(lines)


def character_matrix(
    S_list: list[np.ndarray],
    reference_index: int,
) -> np.ndarray:
    """Return the character of a single reference state across the path.

    For each step, return the |S[reference_index, j]|^2 for each new state j.
    Shape: (n_steps, n_new).

    This is what the paper's Figure 4a plots (state character vs geometry).
    """
    if not S_list:
        return np.zeros((0, 0))
    n_new = S_list[0].shape[1]
    out = np.zeros((len(S_list), n_new))
    for step, S in enumerate(S_list):
        out[step, :] = np.abs(S[reference_index, :]) ** 2
    return out


# ---------------------------------------------------------------- plotting

def plot_character(
    char: np.ndarray,
    path_coords: np.ndarray,
    reference_label: str = "reference state",
    outpath: str | None = None,
):
    """Plot the character of a reference state across a displacement path.

    Parameters
    ----------
    char : (n_steps, n_new) array of |<ref|S_j>|^2
    path_coords : (n_steps,) coordinate along the path (e.g. bond length)
    reference_label : label for the plot title
    outpath : if not None, save the figure to this path
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 4))
    for j in range(char.shape[1]):
        ax.plot(path_coords, char[:, j], marker="o", label=f"S{j+1}")
    ax.set_xlabel("displacement coordinate")
    ax.set_ylabel(f"|<{reference_label}|S_j>|²")
    ax.set_title(f"Character of {reference_label} across the path")
    ax.set_ylim(-0.02, 1.02)
    ax.legend(loc="best", fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    if outpath:
        fig.savefig(outpath, dpi=150)
    plt.close(fig)
