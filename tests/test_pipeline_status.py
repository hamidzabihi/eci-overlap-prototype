"""
Tests for the pipeline status reporter.
"""

from src.pipeline_status import (
    pipeline_stages,
    missing_stages,
    implemented_stages,
    format_pipeline_report,
)


def test_stages_are_numbered_sequentially():
    stages = pipeline_stages()
    nums = [s.num for s in stages]
    assert nums == sorted(nums)
    assert nums == list(range(1, len(stages) + 1))


def test_statuses_are_valid():
    valid = {"present", "done", "partial", "missing"}
    for s in pipeline_stages():
        assert s.status in valid, f"{s.name} has invalid status {s.status}"


def test_expected_missing_stages():
    """After Phase 5.3, stages 6 and 7 are done. Stages 8..11 remain missing."""
    missing = [s.num for s in missing_stages()]
    assert missing == [8, 9, 10, 11]


def test_expected_implemented_stages():
    """Stages 1..7 are present or done after Phase 5.3."""
    implemented = [s.num for s in implemented_stages()]
    assert implemented == [1, 2, 3, 4, 5, 6, 7]


def test_report_mentions_deliverables():
    report = format_pipeline_report()
    assert "wavefunction overlaps" in report
    assert "energy gradients" in report
    assert "spin-orbit couplings" in report
    assert "Missing sections" not in report  # this is the pipeline report, not the QM.out gap map


def test_report_has_summary_counts():
    report = format_pipeline_report()
    assert "implemented (present or done):  7" in report
    assert "missing:                        4" in report
