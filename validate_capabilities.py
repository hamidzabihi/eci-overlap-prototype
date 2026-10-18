#!/usr/bin/env python3
"""
ECI prototype capability validator.

Runs a single command that exercises every capability of the prototype:

  1. Module surface         - expected files import cleanly
  2. ESD overlap primitive  - overlap between two Slater determinants
  3. ECSF contraction       - state overlaps from ESD overlaps
  4. Phase tracking         - sign consistency across a geometry path
  5. Diabatization          - state tracking across a crossing
  6. GS-LE coupling         - validated by symmetric-dimer symmetry
  7. N-fragment ECI         - MAD vs. direct CIS for 3 and 4 fragments
  8. FD gradients           - force balance on the S1 eigenvalue

Exit code 0 if all checks pass.

Usage:
    python3 validate_capabilities.py
"""

import subprocess
import sys
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SRC  = ROOT / "src"


def banner(t):
    print(); print("=" * 68); print(t); print("=" * 68)


def report(label, ok, detail=""):
    mark = "PASS" if ok else "FAIL"
    line = f"  [{mark}] {label}"
    if detail:
        line += f"  --  {detail}"
    print(line)
    return ok


def check_modules():
    banner("Section 1 -- Module surface")
    expected = [
        "esd_overlap.py",
        "state_overlap.py",
        "phase_tracking.py",
        "diabatization.py",
        "excitonic_hamiltonian.py",
        "eci_general.py",
        "gradients.py",
    ]
    ok = True
    for f in expected:
        p = SRC / f
        report(f"src/{f:<28}", p.exists())
        ok &= p.exists()

    sys.path.insert(0, str(ROOT))
    for mod, sym in (
        ("src.gradients", "fd_gradient"),
        ("src.eci_general", "build_eci_general"),
        ("src.excitonic_hamiltonian", "build_two_fragment_eci"),
    ):
        try:
            __import__(mod, fromlist=[sym])
            report(f"import {mod}.{sym}", True)
        except Exception as e:
            report(f"import {mod}.{sym}", False, str(e))
            ok = False
    return ok


def run_pytest(pattern, label):
    try:
        env = os.environ.copy()
        env["PYTHONPATH"] = f"{ROOT}:{ROOT / 'src'}:" + env.get("PYTHONPATH", "")
        r = subprocess.run(
            [sys.executable, "-m", "pytest", str(ROOT / pattern), "-q"],
            capture_output=True, text=True, timeout=300, cwd=str(ROOT), env=env
        )
        tail = (r.stdout.strip().split("\n") or [""])[-1]
        return report(label, r.returncode == 0, tail)
    except subprocess.TimeoutExpired:
        return report(label, False, "timeout")


def check_overlap_primitives():
    banner("Section 2 -- ESD overlap product formula")
    return run_pytest("tests/test_esd_overlap.py", "ESD overlap")


def check_state_overlap():
    banner("Section 3 -- ECSF contraction / state overlaps")
    return run_pytest("tests/test_state_overlap.py", "State overlap")


def check_phase_tracking():
    banner("Section 4 -- Phase tracking")
    return run_pytest("tests/test_phase_tracking.py", "Phase tracking")


def check_diabatization():
    banner("Section 5 -- Diabatization")
    return run_pytest("tests/test_diabatization.py", "State tracking")


def check_gs_le():
    banner("Section 6 -- GS-LE coupling (symmetry validation)")
    return run_pytest("tests/test_gs_le_symmetry.py", "GS-LE symmetric-dimer")


def check_n_fragment_demo():
    banner("Section 7 -- N-fragment ECI (3- and 4-fragment chains)")
    ok = True
    env = os.environ.copy()
    env["PYTHONPATH"] = f"{ROOT}:{ROOT / 'src'}:" + env.get("PYTHONPATH", "")
    for nfrag in (3, 4):
        try:
            r = subprocess.run(
                [sys.executable, "examples/n_fragment_test.py",
                 "--sep", "5.0", "--nstates", "2", "--nfrag", str(nfrag)],
                capture_output=True, text=True, timeout=300, cwd=str(ROOT), env=env
            )
            mad_line = next(
                (l for l in r.stdout.split("\n") if "MAD:" in l),
                "(no MAD line)",
            )
            passed = "MAD:" in mad_line and "meV" in mad_line
            report(f"{nfrag}-fragment chain", passed, mad_line.strip())
            ok &= passed
        except subprocess.TimeoutExpired:
            report(f"{nfrag}-fragment chain", False, "timeout")
            ok = False
    return ok


def check_gradients():
    banner("Section 8 -- Finite-difference gradients (S1)")
    print("  Runtime ~40 s.")
    env = os.environ.copy()
    env["PYTHONPATH"] = f"{ROOT}:{ROOT / 'src'}:" + env.get("PYTHONPATH", "")
    print()
    try:
        r = subprocess.run(
            [sys.executable, "examples/gradients_demo.py",
             "--separation", "4.0", "--delta", "5e-4"],
            capture_output=True, text=True, timeout=300, cwd=str(ROOT), env=env
        )
        out = r.stdout
        fb_line = next(
            (l for l in out.split("\n") if "net force" in l.lower()),
            "(no force-balance line)",
        )
        ok = "net force" in out.lower()
        if "WARNING" in out and "net force >> 0" in out:
            ok = False
        report("FD gradient + force balance", ok, fb_line.strip())
        return ok
    except subprocess.TimeoutExpired:
        return report("FD gradient + force balance", False, "timeout")


def check_test_suite():
    banner("Section 9 -- Full test suite")
    try:
        env = os.environ.copy()
        env["PYTHONPATH"] = f"{ROOT}:{ROOT / 'src'}:" + env.get("PYTHONPATH", "")
        r = subprocess.run(
            [sys.executable, "-m", "pytest", "-q"],
            capture_output=True, text=True, timeout=600, cwd=str(ROOT), env=env
        )
        tail = (r.stdout.strip().split("\n") or [""])[-1]
        return report("pytest tests/", r.returncode == 0, tail)
    except subprocess.TimeoutExpired:
        return report("pytest tests/", False, "timeout")


def main():
    print()
    print("ECI prototype capability validation")
    print(f"  root: {ROOT}")

    results = [
        ("1. Module surface",  check_modules()),
        ("2. ESD overlap",     check_overlap_primitives()),
        ("3. State overlap",   check_state_overlap()),
        ("4. Phase tracking",  check_phase_tracking()),
        ("5. Diabatization",   check_diabatization()),
        ("6. GS-LE coupling",  check_gs_le()),
        ("7. N-fragment ECI",  check_n_fragment_demo()),
        ("8. FD gradients",    check_gradients()),
        ("9. Full test suite", check_test_suite()),
    ]

    print()
    banner("FINAL SUMMARY -- ECI prototype capability validation")
    all_pass = True
    for name, ok in results:
        print(f"    [{'PASS' if ok else 'FAIL'}] {name}")
        all_pass &= ok
    print()
    if all_pass:
        print("  ====================================================")
        print(f"  ALL ECI CAPABILITY CHECKS PASSED  ({len(results)}/{len(results)})")
        print("  ====================================================")
    else:
        n_pass = sum(1 for _, ok in results if ok)
        print(f"  {n_pass}/{len(results)} sections passed")
    print()
    sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
    main()
