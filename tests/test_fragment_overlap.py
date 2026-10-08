"""
Tests for the fragment-overlap file writers.

These tests validate the file formats in isolation, without running
wfoverlap. They check headers, dimensions, and the round-trip properties
of the formats.
"""

from pathlib import Path
import numpy as np
import pytest

from src.fragment_overlap import (
    write_mo_file,
    write_dets_file_ci,
    write_ao_overlap_file,
    write_cioverlap_input,
    WfoverlapFiles,
)


# ------------------------------------------------------------ MO file

def test_mo_file_header(tmp_path):
    mo = np.array([[0.5, 0.5], [0.5, -0.5]])
    out = tmp_path / "mos"
    write_mo_file(out, mo)
    lines = out.read_text().splitlines()
    assert lines[0] == "2mocoef"
    assert lines[1].strip() == "header"
    # dimension line: "     2   2"
    dim_line = lines[5].split()
    assert dim_line == ["2", "2"]


def test_mo_file_contains_all_coefficients(tmp_path):
    mo = np.array([[0.5, 0.5, 0.1], [0.5, -0.5, 0.2], [0.3, 0.4, 0.9]])
    out = tmp_path / "mos"
    write_mo_file(out, mo)
    text = out.read_text()
    # All values should appear as " 5.0000000000000000e-01" style
    for x in mo.flatten():
        assert f"{x:.16e}" in text or f" {x:.16e} " in text


def test_mo_file_orbocc_section(tmp_path):
    mo = np.eye(2)
    out = tmp_path / "mos"
    write_mo_file(out, mo)
    text = out.read_text()
    assert "orbocc" in text
    assert "(*)" in text


# ------------------------------------------------------------ dets file

def test_dets_file_header(tmp_path):
    dets = ["dd", "ee"]
    coefs = np.array([[1.0, 0.0], [0.0, 1.0]])
    out = tmp_path / "dets"
    write_dets_file_ci(out, dets, coefs, norb=2)
    header = out.read_text().splitlines()[0]
    assert header == "2 2 2"


def test_dets_file_wrong_length_raises(tmp_path):
    with pytest.raises(ValueError):
        write_dets_file_ci(tmp_path / "dets", ["ddd"], np.array([[1.0]]), norb=2)


def test_dets_file_sorting(tmp_path):
    dets = ["ee", "dd", "aa"]
    coefs = np.eye(3)
    out = tmp_path / "dets"
    write_dets_file_ci(out, dets, coefs, norb=2)
    lines = out.read_text().splitlines()
    # After header, the determinant rows should be in reverse lex order
    # with mapping e=0, a=1, b=2, d=3:
    #   dd -> (3,3)
    #   aa -> (1,1)
    #   ee -> (0,0)
    # reverse order: dd, aa, ee
    row_dets = [ln[:2] for ln in lines[1:]]
    assert row_dets == ["dd", "aa", "ee"]


def test_dets_file_roundtrip(tmp_path):
    """Write and read back; verify header and determinant set."""
    dets = ["dd", "da", "ad", "aa"]
    coefs = np.array([
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 0.7, 0.7, 0.0],
    ])
    out = tmp_path / "dets"
    write_dets_file_ci(out, dets, coefs, norb=2)
    lines = out.read_text().splitlines()
    assert lines[0] == "2 2 4"
    assert len(lines) == 5  # header + 4 rows
    written_dets = {row[:2] for row in lines[1:]}
    assert written_dets == set(dets)


# ------------------------------------------------------------ AO overlap

def test_ao_overlap_header(tmp_path):
    sao = np.eye(3)
    out = tmp_path / "ao"
    write_ao_overlap_file(out, sao)
    lines = out.read_text().splitlines()
    assert lines[0] == "3 3"
    assert len(lines) == 4


def test_ao_overlap_non_square_raises(tmp_path):
    with pytest.raises(ValueError):
        write_ao_overlap_file(tmp_path / "ao", np.ones((2, 3)))


def test_ao_overlap_roundtrip(tmp_path):
    rng = np.random.default_rng(42)
    sao = rng.random((4, 4))
    out = tmp_path / "ao"
    write_ao_overlap_file(out, sao)
    lines = out.read_text().splitlines()
    matrix = np.array([[float(x) for x in line.split()] for line in lines[1:]])
    np.testing.assert_allclose(matrix, sao, atol=1e-15)


# ------------------------------------------------------------ cioverlap.input

def test_cioverlap_input_format(tmp_path):
    files = WfoverlapFiles(
        mo_a="mos_a", mo_b="mos_b",
        det_a="dets_a", det_b="dets_b",
        ao_overlap="AO_overl.mixed",
    )
    out = tmp_path / "cioverlap.input"
    write_cioverlap_input(out, files)
    lines = out.read_text().splitlines()
    assert "a_mo=mos_a" in lines
    assert "b_mo=mos_b" in lines
    assert "a_det=dets_a" in lines
    assert "b_det=dets_b" in lines
    assert "mix_aoovl=AO_overl.mixed" in lines
    assert "ao_read=0" in lines
