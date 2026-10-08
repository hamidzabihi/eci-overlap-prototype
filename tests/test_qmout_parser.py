"""
Validate the QM.out parser against the SHARC_ECI example output.
"""

from pathlib import Path
import numpy as np
import pytest

from src.qmout_parser import (
    parse_qmout, gap_report, SHARC_NAMD_REQUIRED,
)

ECI_QMout = Path(__file__).parent / "fixtures" / "eci_bodipy_dimer.out"


@pytest.fixture(scope="module")
def eci_data():
    if not ECI_QMout.is_file():
        pytest.skip("ECI QM.out not found at " + str(ECI_QMout))
    return parse_qmout(ECI_QMout)


def test_basic_info(eci_data):
    assert eci_data.states == [6, 0, 6]
    assert eci_data.charges == [0, 1, 0]
    assert eci_data.nmstates == 24
    assert eci_data.natom == 42
    assert eci_data.npc == 0


def test_hamiltonian_shape_and_symmetry(eci_data):
    assert eci_data.h is not None
    assert eci_data.h.shape == (24, 24)
    np.testing.assert_allclose(eci_data.h, eci_data.h.conj().T, atol=1e-10)


def test_dipole_matrices_shape(eci_data):
    assert eci_data.dm is not None
    assert eci_data.dm.shape == (3, 24, 24)
    assert abs(eci_data.dm[0, 0, 3] - 1.061984309357) < 1e-8


def test_runtime_parsed(eci_data):
    assert eci_data.runtime is not None
    assert abs(eci_data.runtime - 233.0) < 1e-6


def test_gap_map_shows_missing_namd_sections(eci_data):
    missing = SHARC_NAMD_REQUIRED - eci_data.present_sections
    assert 3 in missing
    assert 5 in missing
    assert 6 in missing
    assert 7 in missing
    assert 13 in missing


def test_gap_report_is_readable(eci_data):
    report = gap_report(eci_data)
    assert "Missing sections for SHARC NAMD" in report
    assert "Hamiltonian Matrix" in report
    assert "Gradient Vectors" in report
