#!/usr/bin/env python3
"""Check sign-reflection and axis-exchange twist symmetry with exact ED."""

from __future__ import annotations

import argparse
import csv
import subprocess
import tempfile
from pathlib import Path


def run(ed: Path, output: Path, phix: float, phiy: float):
    result = subprocess.run(
        [
            str(ed),
            "--lx", "2", "--ly", "2",
            "--u", "7", "--tx", "1", "--ty", "1", "--tp", "0", "--v", "0",
            "--phix", str(phix), "--phiy", str(phiy),
            "--beta", "12.5", "--mu-min", "-3", "--mu-max", "4", "--mu-count", "41",
            "--output", str(output),
        ],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if result.returncode:
        raise AssertionError(result.stderr or result.stdout)


def rows(path: Path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def compare(reference: Path, candidate: Path, tolerance: float = 2e-10):
    left = rows(reference)
    right = rows(candidate)
    assert len(left) == len(right)
    columns = ("mu", "n", "charge_correlation", "compressibility", "log_partition")
    for index, (a, b) in enumerate(zip(left, right)):
        for column in columns:
            error = abs(float(a[column]) - float(b[column]))
            assert error < tolerance, (index, column, error)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ed", required=True, type=Path)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="ftlm-twist-d4-") as directory:
        root = Path(directory)
        cases = {
            "base": (0.1875, 0.3125),
            "negx": (-0.1875, 0.3125),
            "negy": (0.1875, -0.3125),
            "swap": (0.3125, 0.1875),
            "boundary_pos": (0.5, 0.15),
            "boundary_neg": (-0.5, 0.15),
        }
        for name, (phix, phiy) in cases.items():
            run(args.ed, root / f"{name}.csv", phix, phiy)
        for name in ("negx", "negy", "swap"):
            compare(root / "base.csv", root / f"{name}.csv")
        compare(root / "boundary_pos.csv", root / "boundary_neg.csv")
    print("PASS ED twist symmetry: signs, kx<->ky, and periodic pi/4 boundary")


if __name__ == "__main__":
    main()
