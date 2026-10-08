from pathlib import Path
import shutil
import numpy as np
import pytest

from src.run_wfoverlap import (
    run_wfoverlap, parse_overlap_matrix, parse_reference_output
)

SHARC_ROOT = Path.home() / "sharc4"
WATER = SHARC_ROOT / "wfoverlap" / "test_jobs" / "water"
WFOVERLAP = SHARC_ROOT / "bin" / "wfoverlap.x"
REF_OUT = WATER / "REF_FILES" / "ciovl.out.1thr"
IN_FILES = WATER / "IN_FILES"


@pytest.fixture
def water_workdir(tmp_path):
    workdir = tmp_path / "water"
    shutil.copytree(IN_FILES, workdir)
    return workdir


def test_wfoverlap_binary_exists():
    assert WFOVERLAP.is_file(), f"wfoverlap.x not found at {WFOVERLAP}"
    assert REF_OUT.is_file(), f"reference not found at {REF_OUT}"


def test_water_overlap_matches_reference(water_workdir):
    stdout = run_wfoverlap("ciovl.in", executable=WFOVERLAP, cwd=water_workdir)
    S = parse_overlap_matrix(stdout, renormalized=False)
    S_ref = parse_reference_output(REF_OUT, renormalized=False)
    np.testing.assert_allclose(S, S_ref, atol=1e-8, rtol=0)


def test_water_renormalized_overlap_matches_reference(water_workdir):
    stdout = run_wfoverlap("ciovl.in", executable=WFOVERLAP, cwd=water_workdir)
    S = parse_overlap_matrix(stdout, renormalized=True)
    S_ref = parse_reference_output(REF_OUT, renormalized=True)
    np.testing.assert_allclose(S, S_ref, atol=1e-8, rtol=0)


def test_water_overlap_is_physical(water_workdir):
    stdout = run_wfoverlap("ciovl.in", executable=WFOVERLAP, cwd=water_workdir)
    S = parse_overlap_matrix(stdout, renormalized=False)
    assert np.all(np.abs(S) <= 1.0 + 1e-10)
