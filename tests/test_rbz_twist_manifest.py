#!/usr/bin/env python3
"""Validate the pi/40 Gamma-centered rBZ grid and D4 representatives."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


GRID_SIZE = 20


def orbit(ix: int, iy: int) -> set[tuple[int, int]]:
    points = set()
    for sx in (-1, 1):
        for sy in (-1, 1):
            points.add(((sx * ix) % GRID_SIZE, (sy * iy) % GRID_SIZE))
            points.add(((sx * iy) % GRID_SIZE, (sy * ix) % GRID_SIZE))
    return points


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--campaign-dir", required=True, type=Path)
    args = parser.parse_args()
    directory = args.campaign_dir.resolve()
    with (directory / "twists.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    metadata = json.loads((directory / "twist_quadrature.json").read_text())

    assert metadata["tag"] == "rbz20_gamma_dk_pi40"
    assert metadata["cluster_lx"] == metadata["cluster_ly"] == 4
    assert metadata["grid_nx"] == metadata["grid_ny"] == GRID_SIZE
    assert metadata["representatives"] == len(rows) == 66
    assert metadata["effective_twists"] == 400
    assert metadata["includes_zero"] is True
    assert math.isclose(metadata["spacing_over_pi"], 1 / 40)

    expected_pairs = {(ix, iy) for iy in range(11) for ix in range(iy + 1)}
    observed_pairs = set()
    expanded = set()
    total_weight = 0
    next_twist_id = 100
    for row in rows:
        ix = int(row["grid_ix"])
        iy = int(row["grid_iy"])
        if (ix, iy) == (0, 0):
            twist_id = 0
            expected_seed = 12345
        elif (ix, iy) == (5, 5):
            twist_id = 5
            expected_seed = 5012360
        else:
            twist_id = next_twist_id
            expected_seed = 12345 + 1000003 * twist_id
            next_twist_id += 1
        assert row["twist_id"] == f"{twist_id:03d}"
        assert int(row["seed"]) == expected_seed
        x = float(row["phix"])
        y = float(row["phiy"])
        assert 0 <= ix <= iy <= 10
        assert math.isclose(x, ix / GRID_SIZE, abs_tol=1e-12)
        assert math.isclose(y, iy / GRID_SIZE, abs_tol=1e-12)
        observed_pairs.add((ix, iy))
        assert math.isclose(float(row["kx_over_pi"]), ix / 40, abs_tol=1e-12)
        assert math.isclose(float(row["ky_over_pi"]), iy / 40, abs_tol=1e-12)
        points = orbit(ix, iy)
        expected_weight = len(points)
        assert int(row["multiplicity"]) == expected_weight
        total_weight += expected_weight
        assert expanded.isdisjoint(points)
        expanded.update(points)

    assert observed_pairs == expected_pairs
    assert next_twist_id == 164
    assert total_weight == len(expanded) == 400
    assert expanded == {(ix, iy) for ix in range(20) for iy in range(20)}
    special = {(int(row["grid_ix"]), int(row["grid_iy"])): int(row["multiplicity"])
               for row in rows}
    assert special[(0, 0)] == special[(10, 10)] == 1
    assert special[(0, 10)] == 2
    assert special[(0, 1)] == special[(1, 1)] == special[(1, 10)] == 4
    assert special[(1, 2)] == 8
    print("PASS representatives=66 effective_twists=400 spacing=pi/40 zero=included")


if __name__ == "__main__":
    main()
