"""
Compute fragment-level wavefunction overlaps via the SHARC wfoverlap engine.

This module produces the three input files that wfoverlap reads:

    mos         - MO coefficients      (format: 2mocoef, from SHARC_GAUSSIAN)
    dets        - determinants + CI    (format: nstates norb ndets + rows)
    AO_overl.mixed - mixed AO overlap  (format: nAO nAO + matrix)

and orchestrates the wfoverlap call. The file formats are documented
in src/README_formats.md and match those produced by SHARC_GAUSSIAN.

References
----------
- sharc4/bin/SHARC_GAUSSIAN.py: get_MO_from_chk (MO format)
- sharc4/bin/SHARC_GAUSSIAN.py: format_ci_vectors (dets format)
- sharc4/bin/SHARC_GAUSSIAN.py: _create_aoovl (AO overlap format)
- sharc4/wfoverlap/source/inputmod.f90 (cioverlap.input keys)
"""

from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence
import numpy as np


# ------------------------------------------------------------------ MO file

MO_COEF_FMT = " % 12.16e "
MO_OCC_FMT = " % 6.12e "


def write_mo_file(path: Path | str, mo_coeff: np.ndarray, unrestricted: bool = False, mo_coeff_b: np.ndarray | None = None) -> None:
    """Write a MO coefficient file in the SHARC 2mocoef format.

    Parameters
    ----------
    path         : output path
    mo_coeff     : (nAO, nMO) array of alpha MO coefficients (or both for restricted)
    unrestricted : if True, also write beta MOs from mo_coeff_b
    mo_coeff_b   : (nAO, nMO_beta) array of beta MO coefficients

    Notes
    -----
    The format is copied from SHARC_GAUSSIAN.get_MO_from_chk: header, MO block
    (3 coefficients per line, format " % 12.16e "), then orbocc block.
    """
    path = Path(path)
    nao, nmo_a = mo_coeff.shape
    if unrestricted and mo_coeff_b is not None:
        nmo = nmo_a + mo_coeff_b.shape[1]
    else:
        nmo = nmo_a

    lines = []
    lines.append("2mocoef")
    lines.append("    header")
    lines.append("     1")
    lines.append("    MO-coefficients from Python")
    lines.append("     1")
    lines.append(f"     {nao}   {nmo}")
    lines.append("     a")
    lines.append("    mocoef")
    lines.append("    (*)")

    # Alpha MO coefficients
    x = 0
    buf = ""
    for imo in range(mo_coeff.shape[1]):
        for iao in range(nao):
            if x >= 3:
                lines.append(buf)
                buf = ""
                x = 0
            buf += MO_COEF_FMT % mo_coeff[iao, imo]
            x += 1
        if x > 0:
            lines.append(buf)
            buf = ""
            x = 0

    # Beta MO coefficients (if unrestricted)
    if unrestricted and mo_coeff_b is not None:
        for imo in range(mo_coeff_b.shape[1]):
            for iao in range(nao):
                if x >= 3:
                    lines.append(buf)
                    buf = ""
                    x = 0
                buf += MO_COEF_FMT % mo_coeff_b[iao, imo]
                x += 1
            if x > 0:
                lines.append(buf)
                buf = ""
                x = 0

    # Orbital occupations (all zero, as in SHARC_GAUSSIAN)
    lines.append("orbocc")
    lines.append("(*)")
    x = 0
    buf = ""
    for i in range(nmo):
        if x >= 3:
            lines.append(buf)
            buf = ""
            x = 0
        buf += MO_OCC_FMT % 0.0
        x += 1
    if x > 0:
        lines.append(buf)

    path.write_text("\n".join(lines) + "\n")


# ------------------------------------------------------------------ dets file

DETS_COEF_FMT = " %15.10E "

_SYMBOL_TO_INT = {"e": 0, "a": 1, "b": 2, "d": 3}


def write_dets_file_ci(
    path: Path | str,
    detstrings: Sequence[str],
    coefficients: np.ndarray,
    norb: int,
) -> None:
    """Write a dets file in the format wfoverlap reads.

    Parameters
    ----------
    path         : output path
    detstrings   : list of detstrings, each of length norb
    coefficients : (nstate, nsd) array; coefficient[s, d] is the coefficient
                   of determinant d in state s
    norb         : number of orbitals (length of each detstring)

    Notes
    -----
    Determinants are sorted in reverse lexicographic order, matching
    SHARC_GAUSSIAN.format_ci_vectors. Determinants are written first,
    then nstates coefficients, formatted as " %15.10E ".
    """
    path = Path(path)
    if coefficients.ndim != 2:
        raise ValueError("coefficients must be 2D (nstate, nsd)")
    nstate, nsd = coefficients.shape
    if nsd != len(detstrings):
        raise ValueError(f"coefficients has {nsd} columns, detstrings has {len(detstrings)}")
    for i, d in enumerate(detstrings):
        if len(d) != norb:
            raise ValueError(f"detstring {i} has length {len(d)}, expected {norb}")

    # Sort determinants by their integer-tuple representation, reverse
    order = sorted(
        range(nsd),
        key=lambda i: tuple(_SYMBOL_TO_INT[c] for c in detstrings[i]),
        reverse=True,
    )

    lines = [f"{nstate} {norb} {nsd}"]
    for d in order:
        row = detstrings[d]
        for s in range(nstate):
            row += DETS_COEF_FMT % coefficients[s, d]
        lines.append(row)
    path.write_text("\n".join(lines) + "\n")


# ------------------------------------------------------------------ AO overlap

def write_ao_overlap_file(path: Path | str, sao_mixed: np.ndarray) -> None:
    """Write the mixed AO overlap matrix in native ascii format.

    Format:
        nAO nAO
        row1 (nAO values)
        ...
        row_nAO

    Parameters
    ----------
    path      : output path
    sao_mixed : (nAO, nAO) array; block [old, new] of the full overlap matrix
    """
    path = Path(path)
    nao = sao_mixed.shape[0]
    if sao_mixed.shape[1] != nao:
        raise ValueError(f"AO overlap must be square, got {sao_mixed.shape}")
    lines = [f"{nao} {nao}"]
    for row in sao_mixed:
        lines.append(" ".join(f"{x: .15e}" for x in row))
    path.write_text("\n".join(lines) + "\n")


# ------------------------------------------------------------------ wfoverlap input

@dataclass
class WfoverlapFiles:
    """Paths to the three files wfoverlap reads."""
    mo_a: str
    mo_b: str
    det_a: str
    det_b: str
    ao_overlap: str


def write_cioverlap_input(path: Path | str, files: WfoverlapFiles) -> None:
    """Write a cioverlap.input file pointing at the given data files.

    Format (from inputmod.f90):
        a_mo=...        MO coefficients at geometry A
        b_mo=...        MO coefficients at geometry B
        a_det=...       determinants at A
        b_det=...       determinants at B
        mix_aoovl=...   mixed AO overlap
        ao_read=0       native ascii format
    """
    path = Path(path)
    lines = [
        f"a_mo={files.mo_a}",
        f"b_mo={files.mo_b}",
        f"a_det={files.det_a}",
        f"b_det={files.det_b}",
        f"mix_aoovl={files.ao_overlap}",
        "ao_read=0",
    ]
    path.write_text("\n".join(lines) + "\n")


# ------------------------------------------------------------------ orchestration

def compute_fragment_overlap(
    mo_a: str,
    mo_b: str,
    det_a: str,
    det_b: str,
    ao_overlap: str,
    workdir: Path | str,
    wfoverlap_exe: Path | str,
    keep_inputs: bool = False,
) -> np.ndarray:
    """Run wfoverlap on the given files and return the overlap matrix.

    All paths are relative to `workdir`. This function writes the
    cioverlap.input file and executes wfoverlap.x.

    Parameters
    ----------
    mo_a, mo_b, det_a, det_b, ao_overlap : filenames in workdir
    workdir       : directory containing the input files
    wfoverlap_exe : path to wfoverlap.x
    keep_inputs   : if False, cioverlap.input is deleted after the run

    Returns
    -------
    (nstate_A, nstate_B) complex overlap matrix
    """
    from src.run_wfoverlap import run_wfoverlap, parse_overlap_matrix

    workdir = Path(workdir)
    input_file = workdir / "cioverlap.input"
    files = WfoverlapFiles(
        mo_a=mo_a, mo_b=mo_b, det_a=det_a, det_b=det_b, ao_overlap=ao_overlap,
    )
    write_cioverlap_input(input_file, files)
    try:
        stdout = run_wfoverlap("cioverlap.input", executable=wfoverlap_exe, cwd=workdir)
    finally:
        if not keep_inputs and input_file.exists():
            input_file.unlink()
    return parse_overlap_matrix(stdout, renormalized=False)
