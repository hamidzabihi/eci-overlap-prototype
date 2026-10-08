"""
Lightweight parser for SHARC QM.out files.

Only parses the sections that the ECI interface currently produces:
    ! 0  Basic information
    ! 1  Hamiltonian Matrix
    ! 2  Dipole Moment Matrices
    ! 8  Runtime

The full SHARC spec has many more sections (see SECTIONS below).
This module intentionally does NOT use SHARC's own QMout class, so that
the prototype stays self-contained and can parse partial QM.out files.
"""

from __future__ import annotations
import re
from pathlib import Path
from dataclasses import dataclass, field
import numpy as np


# ---------------------------------------------------------------- spec

# Full SHARC QM.out section specification, taken from
# sharc4/lib/qmout.py (QMout.write / writeQMoutXXX methods).
SECTIONS = {
    0:   "Basic information",
    1:   "Hamiltonian Matrix",
    2:   "Dipole Moment Matrices",
    3:   "Gradient Vectors",
    5:   "Non-adiabatic couplings (ddr)",
    6:   "Overlap matrix",
    7:   "Wave function phases",
    8:   "Runtime",
    12:  "Dipole moment derivatives",
    13:  "Spin-Orbit coupling derivatives",
    20:  "Property Matrices",
    21:  "Property Vectors",
    22:  "Atomwise multipolar density fits",
    23:  "Property Scalars",
    24:  "Total/Spin/Partial 1-particle density matrices",
    25:  "Mole PySCF object",
    30:  "Point Charge Gradient Vectors",
    31:  "Non-adiabatic couplings on point charges (ddr)",
    32:  "Dipole moment derivatives on point charges",
    33:  "Spin-Orbit coupling derivatives on point charges",
    41:  "Magnetic Dipole Moment Matrices",
    42:  "Electric Quadrupole Moment Matrices",
    999: "Notes",
}

# Sections required for surface-hopping NAMD in SHARC.
# Derived from QMout.allocate(): h, grad, nacdr/overlap, phases, socdr.
SHARC_NAMD_REQUIRED = {1, 3, 5, 6, 7, 13}


# ---------------------------------------------------------------- parser

@dataclass
class QMoutData:
    """Container for parsed QM.out sections."""
    path: Path
    states: list = field(default_factory=list)
    charges: list = field(default_factory=list)
    nmstates: int = 0
    natom: int = 0
    npc: int = 0
    nstates: int = 0
    h: object = None
    dm: object = None
    runtime: object = None
    present_sections: set = field(default_factory=set)
    section_headers: dict = field(default_factory=dict)


def _parse_section_header(line):
    """Parse '! <n> <title>' into (n, title)."""
    if not line.startswith("!"):
        return None
    parts = line[1:].strip().split(maxsplit=2)
    if len(parts) < 1:
        return None
    try:
        n = int(parts[0])
    except ValueError:
        return None
    title = parts[1] if len(parts) > 1 else ""
    title = re.sub(r"\s*\(.*\)\s*$", "", title)
    return n, title


def _read_matrix_complex(lines, dim, skip=1):
    """Parse a dim x dim complex matrix, real/imag interleaved on each row.

    The QM.out format writes a dimension line (e.g. "24 24") before the
    matrix data, so by default the first `skip` lines of `lines` are ignored.
    """
    out = np.zeros((dim, dim), dtype=complex)
    for i in range(dim):
        toks = lines[skip + i].split()
        if len(toks) < 2 * dim:
            raise ValueError(
                "row %d has %d tokens, expected %d" % (i, len(toks), 2 * dim)
            )
        for j in range(dim):
            re_ = float(toks[2 * j])
            im_ = float(toks[2 * j + 1])
            out[i, j] = complex(re_, im_)
    return out


def parse_qmout(path):
    """Parse a SHARC QM.out file. Only extracts basic info, Hamiltonian,
    dipole matrices, and runtime. Missing sections are simply not parsed."""
    path = Path(path)
    text = path.read_text()
    lines = text.splitlines()

    data = QMoutData(path=path)

    i = 0
    n = len(lines)
    while i < n:
        header = _parse_section_header(lines[i])
        if header is None:
            i += 1
            continue
        sec, title = header
        data.present_sections.add(sec)
        data.section_headers[sec] = title

        j = i + 1
        block = []
        while j < n and not lines[j].startswith("!"):
            block.append(lines[j])
            j += 1

        while block and block[-1].strip() == "":
            block.pop()

        if sec == 0:
            for line in block:
                line = line.strip()
                if not line:
                    continue
                k, v = line.split(maxsplit=1)
                if k == "states":
                    data.states = [int(x) for x in v.split()]
                elif k == "charges":
                    data.charges = [int(x) for x in v.split()]
                elif k == "nmstates":
                    data.nmstates = int(v)
                elif k == "natom":
                    data.natom = int(v)
                elif k == "npc":
                    data.npc = int(v)
            data.nstates = sum(data.states)

        elif sec == 1:
            if data.nmstates == 0:
                raise ValueError("Hamiltonian section found before basic info")
            dim = data.nmstates
            data.h = _read_matrix_complex(block, dim)

        elif sec == 2:
            if data.nmstates == 0:
                raise ValueError("Dipole section found before basic info")
            dim = data.nmstates
            nblocks = len(block) // (dim + 1)
            data.dm = np.zeros((3, dim, dim), dtype=complex)
            for b in range(nblocks):
                start = b * (dim + 1)
                # each block is: one "dim dim pol X" header + dim data rows
                matrix = _read_matrix_complex(block[start:], dim, skip=1)
                if b < 3:
                    data.dm[b] = matrix

        elif sec == 8:
            if block:
                data.runtime = float(block[0].split()[0])

        i = j

    return data


def gap_report(data):
    """Produce a human-readable report comparing present vs required sections."""
    lines = []
    lines.append("QM.out file: " + str(data.path))
    lines.append("  states:   " + str(data.states))
    lines.append("  charges:  " + str(data.charges))
    lines.append("  nmstates: " + str(data.nmstates))
    lines.append("  natom:    " + str(data.natom))
    lines.append("  npc:      " + str(data.npc))
    lines.append("")
    lines.append("%4s  %-48s %9s %10s" % ("sec", "title", "present", "required"))
    lines.append("-" * 76)
    for sec in sorted(SECTIONS):
        title = SECTIONS[sec]
        present = "yes" if sec in data.present_sections else "-"
        required = "yes" if sec in SHARC_NAMD_REQUIRED else "-"
        lines.append("%4d  %-48s %9s %10s" % (sec, title, present, required))
    lines.append("")
    missing = SHARC_NAMD_REQUIRED - data.present_sections
    if missing:
        lines.append("Missing sections for SHARC NAMD: " + str(sorted(missing)))
    else:
        lines.append("All NAMD-required sections present.")
    return "\n".join(lines)
