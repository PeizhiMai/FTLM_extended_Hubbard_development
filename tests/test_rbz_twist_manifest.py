#!/usr/bin/env python3
"""Validate the 8x8 midpoint rBZ grid and its D4 representatives."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def orbit(x: float, y: float) -> set[tuple[float, float]]:
    points = set()
    for sx in (-1.0, 1.0):
        for sy in (-1.0, 1.0):
            points.add((sx * x, sy * y))
            points.add((sx * y, sy * x))
    return points


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--campaign-dir", required=True, type=Path)
    args = parser.parse_args()
    directory = args.campaign_dir.resolve()
    with (directory / "twists.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    metadata = json.loads((directory / "twist_quadrature.json").read_text())

    assert metadata["tag"] == "rbz8_midpoint"
    assert metadata["cluster_lx"] == metadata["cluster_ly"] == 4
    assert metadata["grid_nx"] == metadata["grid_ny"] == 8
    assert metadata["representatives"] == len(rows) == 10
    assert metadata["effective_twists"] == 64

    positive = {1 / 16, 3 / 16, 5 / 16, 7 / 16}
    expected_pairs = {(x, y) for x in positive for y in positive if x <= y}
    observed_pairs = set()
    expanded = set()
    total_weight = 0
    for offset, row in enumerate(rows):
        twist_id = 800 + offset
        assert row["twist_id"] == f"{twist_id:03d}"
        assert int(row["seed"]) == 12345 + 1000003 * twist_id
        x = float(row["phix"])
        y = float(row["phiy"])
        assert 0.0 < x <= y < 0.5
        observed_pairs.add((x, y))
        assert math.isclose(float(row["kx_over_pi"]), x / 2, abs_tol=1e-12)
        assert math.isclose(float(row["ky_over_pi"]), y / 2, abs_tol=1e-12)
        points = orbit(x, y)
        expected_weight = 4 if math.isclose(x, y) else 8
        assert len(points) == expected_weight
        assert int(row["multiplicity"]) == expected_weight
        total_weight += expected_weight
        assert expanded.isdisjoint(points)
        expanded.update(points)

    assert observed_pairs == expected_pairs
    assert total_weight == len(expanded) == 64
    full_levels = {-7 / 16, -5 / 16, -3 / 16, -1 / 16,
                   1 / 16, 3 / 16, 5 / 16, 7 / 16}
    assert expanded == {(x, y) for x in full_levels for y in full_levels}
    print("PASS representatives=10 effective_twists=64 domain=0<kx<=ky<pi/4")


if __name__ == "__main__":
    main()
