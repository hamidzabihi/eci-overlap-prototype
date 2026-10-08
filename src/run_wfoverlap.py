"""
Wrapper around the SHARC wfoverlap Fortran binary.
"""

from __future__ import annotations
import subprocess
from pathlib import Path
import numpy as np


def run_wfoverlap(input_file, executable, cwd=None, maxmem_mb=None, timeout=600):
    cmd = [str(executable), "-f", str(input_file)]
    if maxmem_mb is not None:
        cmd += ["-m", str(maxmem_mb)]
    result = subprocess.run(cmd, cwd=cwd, capture_output=True,
                            text=True, check=False, timeout=timeout)
    if result.returncode != 0:
        raise RuntimeError(
            f"wfoverlap failed (exit {result.returncode})\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    return result.stdout


def parse_overlap_matrix(stdout, renormalized=False):
    lines = stdout.splitlines()
    target = ("Renormalized overlap matrix" if renormalized
              else "Overlap matrix <PsiA_i|PsiB_j>")
    start = None
    for i, line in enumerate(lines):
        if target in line:
            start = i
            break
    if start is None:
        raise ValueError(f"Could not find '{target}' in output")
    rows = []
    started = False
    for line in lines[start + 2:]:
        s = line.strip()
        if s.startswith("<PsiA"):
            started = True
            body = s.split("|", 1)[1]
            rows.append([float(x) for x in body.split()])
        elif started and s:
            # First <PsiA block has ended. wfoverlap prints three such blocks
            # (raw, renormalized, orthonormalized); we only want the first.
            break
    if not rows:
        raise ValueError("No rows found after header")
    return np.array(rows)


def parse_reference_output(path, renormalized=False):
    return parse_overlap_matrix(Path(path).read_text(), renormalized=renormalized)
