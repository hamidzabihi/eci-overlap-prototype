#!/usr/bin/env python3
"""
Run the ECI QM.out through the parser and print the interface gap map.

Usage:
    python examples/gap_map_demo.py [path/to/QM.out]
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.qmout_parser import parse_qmout, gap_report


def main(path):
    data = parse_qmout(path)
    print(gap_report(data))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        default = Path(__file__).resolve().parent.parent / "tests" / "fixtures" / "eci_bodipy_dimer.out"
        print("No path given, using default: " + str(default))
        main(str(default))
    else:
        main(sys.argv[1])
