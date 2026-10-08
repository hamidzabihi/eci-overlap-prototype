"""
Report on the ECI -> wfoverlap -> SHARC-NAMD pipeline status.

This module encodes a single, factual picture of what exists and what
is missing in the current SHARC-ECI codebase, derived from:

  * sharc4/lib/ECI.py          (the ECI method)
  * sharc4/bin/SHARC_ECI.py    (the SHARC driver)
  * sharc4/wfoverlap/          (the determinant-overlap engine)
  * sharc4/lib/qmout.py        (the QM.out format spec)
  * sharc4/examples/SHARC_ECI/ (a working example)

Everything the reporter claims is verifiable against those files.
"""

from __future__ import annotations
from dataclasses import dataclass, field


# ------------------------------------------------------------------ model

@dataclass
class Stage:
    num: int
    name: str
    status: str                     # "present" | "done" | "missing" | "partial"
    implemented_by: str
    notes: str = ""
    qmout_section: int | None = None  # if this stage produces a QM.out section


# ------------------------------------------------------------------ stages

def pipeline_stages() -> list[Stage]:
    """Return the ordered list of pipeline stages."""
    return [
        Stage(
            num=1,
            name="ECI ESD definition",
            status="present",
            implemented_by="SHARC (lib/ECI.py)",
            notes="excitonic_slater_determinant class; spin adaptation to ECSFs",
        ),
        Stage(
            num=2,
            name="ESD -> detstring mapping",
            status="done",
            implemented_by="this repo (src/esd_to_dets.py, Phase 3)",
            notes="maps FragmentOccupation objects to wfoverlap detstring symbols (d/a/b/e)",
        ),
        Stage(
            num=3,
            name="dets file writer",
            status="done",
            implemented_by="this repo (src/esd_to_dets.py, Phase 3)",
            notes="write_dets_file produces the exact format wfoverlap reads",
        ),
        Stage(
            num=4,
            name="wfoverlap engine",
            status="present",
            implemented_by="SHARC (wfoverlap/source)",
            notes="Fortran determinant-overlap engine, built here with gfortran",
        ),
        Stage(
            num=5,
            name="Python wrapper for wfoverlap",
            status="done",
            implemented_by="this repo (src/run_wfoverlap.py, Phase 1)",
            notes="4 regression tests against SHARC's own water_molcas reference",
        ),
        Stage(
            num=6,
            name="ECSF contraction (O_det -> S_ij)",
            status="missing",
            implemented_by="—",
            notes="contract determinant-level overlap with spin-adaptation matrix U and ECI coeffs",
            qmout_section=6,
        ),
        Stage(
            num=7,
            name="State overlap matrix S_ij",
            status="missing",
            implemented_by="—",
            notes="full-system wavefunction overlaps needed at each NAMD step",
            qmout_section=6,
        ),
        Stage(
            num=8,
            name="Wavefunction phases",
            status="missing",
            implemented_by="—",
            notes="sign/phase tracking for the state-overlap matrix",
            qmout_section=7,
        ),
        Stage(
            num=9,
            name="ECI energy gradients",
            status="missing",
            implemented_by="—",
            notes="derivatives of the effective Hamiltonian w.r.t. nuclear coordinates",
            qmout_section=3,
        ),
        Stage(
            num=10,
            name="Non-adiabatic couplings",
            status="missing",
            implemented_by="—",
            notes="can be computed from overlaps, but the SHARC path is not wired up",
            qmout_section=5,
        ),
        Stage(
            num=11,
            name="Spin-orbit couplings",
            status="missing",
            implemented_by="—",
            notes="fragment SOC integrals + ECSF transformation",
            qmout_section=13,
        ),
    ]


# ------------------------------------------------------------------ report

_STATUS_SYMBOL = {
    "present": "[=]",
    "done":    "[+]",
    "partial": "[~]",
    "missing": "[ ]",
}


def format_pipeline_report() -> str:
    stages = pipeline_stages()

    lines = []
    lines.append("ECI -> wfoverlap -> SHARC NAMD: Pipeline Status")
    lines.append("=" * 78)
    lines.append("")
    lines.append(f"{'#':>2}  {'stage':<34} {'status':<9} {'implemented by':<32}")
    lines.append("-" * 78)

    for s in stages:
        sym = _STATUS_SYMBOL.get(s.status, "[?]")
        lines.append(f"{s.num:>2}  {s.name:<34} {sym} {s.status:<6} {s.implemented_by:<32}")

    lines.append("")

    present = [s for s in stages if s.status in ("present", "done")]
    missing = [s for s in stages if s.status == "missing"]

    lines.append("Summary")
    lines.append("-" * 78)
    lines.append(f"  implemented (present or done): {len(present):>2}")
    lines.append(f"  missing:                       {len(missing):>2}")
    lines.append("")

    if missing:
        lines.append("Missing stages and the ECI2GAME deliverables they map to:")
        lines.append("")
        by_deliverable = {
            "wavefunction overlaps": [6, 7, 8],
            "energy gradients":      [9],
            "non-adiabatic couplings": [10],
            "spin-orbit couplings":  [11],
        }
        for name, nums in by_deliverable.items():
            missing_nums = [s.num for s in missing if s.num in nums]
            if not missing_nums:
                continue
            lines.append(f"  -> {name}")
            for s in stages:
                if s.num in missing_nums:
                    qs = f"  [QM.out section {s.qmout_section}]" if s.qmout_section else ""
                    lines.append(f"       stage {s.num}: {s.name}{qs}")

    lines.append("")
    lines.append("Notes")
    lines.append("-" * 78)
    lines.append("  - Stages 1-5 exist or are implemented in this prototype.")
    lines.append("  - Stages 6-11 are absent from the current SHARC-ECI code.")
    lines.append("  - The QM.out gap map (src/qmout_parser.py) confirms that")
    lines.append("    sections 3, 5, 6, 7, 13 are not written by the ECI interface.")
    lines.append("  - This maps 1:1 onto the ECI2GAME postdoc deliverables.")
    return "\n".join(lines)


def missing_stages() -> list[Stage]:
    return [s for s in pipeline_stages() if s.status == "missing"]


def implemented_stages() -> list[Stage]:
    return [s for s in pipeline_stages() if s.status in ("present", "done")]
