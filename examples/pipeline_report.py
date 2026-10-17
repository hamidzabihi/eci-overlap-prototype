#!/usr/bin/env python3
"""
Print the ECI -> wfoverlap -> SHARC-NAMD pipeline status report.

Usage:
    python examples/pipeline_report.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.pipeline_status import format_pipeline_report


if __name__ == "__main__":
    print(format_pipeline_report())
