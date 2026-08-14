#!/usr/bin/env python3
"""Gate the v2 3x3 FTLM rerun against the existing exact benchmark."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def read(path: Path):
    with path.open(newline="") as stream:
        return [{key: float(value) for key, value in row.items()} for row in csv.DictReader(stream)]


def interpolate(points, x):
    if x <= points[0][0]:
        return points[0][1]
    if x >= points[-1][0]:
        return points[-1][1]
    low, high = 0, len(points) - 1
    while high - low > 1:
        middle = (low + high) // 2
        if points[middle][0] <= x:
            low = middle
        else:
            high = middle
    x0, y0 = points[low]
    x1, y1 = points[high]
    if x1 == x0:
        return 0.5 * (y0 + y1)
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0)


def metric(ed_rows, ftlm_rows):
    ed_rows = sorted(ed_rows, key=lambda row: row["mu"])
    ftlm_rows = sorted(ftlm_rows, key=lambda row: row["mu"])
    if len(ed_rows) != len(ftlm_rows):
        raise RuntimeError("ED and FTLM mu grids have different lengths")
    selected = [
        (ed, ftlm)
        for ed, ftlm in zip(ed_rows, ftlm_rows)
        if 0.65 <= ed["n"] <= 1.0 and 0.65 <= ftlm["n"] <= 1.0
    ]
    ed_curve = sorted((1.0 - ed["n"], ed["compressibility"]) for ed, _ in selected)
    ftlm_curve = sorted((1.0 - ftlm["n"], ftlm["compressibility"]) for _, ftlm in selected)
    x_min = max(0.0, ed_curve[0][0], ftlm_curve[0][0])
    x_max = min(0.35, ed_curve[-1][0], ftlm_curve[-1][0])
    grid = [x_min + (x_max - x_min) * index / 350.0 for index in range(351)]
    differences = [interpolate(ftlm_curve, x) - interpolate(ed_curve, x) for x in grid]
    return {
        "points_fixed_mu": len(selected),
        "x_min": x_min,
        "x_max": x_max,
        "kappa_rms_fixed_x": math.sqrt(sum(value * value for value in differences) / len(differences)),
        "kappa_mean_abs_fixed_x": sum(abs(value) for value in differences) / len(differences),
        "kappa_max_abs_fixed_x": max(abs(value) for value in differences),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ftlm", type=Path, required=True)
    parser.add_argument("--ed-t008", type=Path, required=True)
    parser.add_argument("--ed-t035", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    ftlm = read(args.ftlm)
    cases = [
        (0.08, 12.5, args.ed_t008, 0.015),
        (0.35, 2.857142857142857, args.ed_t035, 0.005),
    ]
    results = []
    for temperature, beta, ed_path, threshold in cases:
        rows = [row for row in ftlm if abs(row["beta"] - beta) < 1e-10]
        result = metric(read(ed_path), rows)
        result.update(
            {
                "temperature": temperature,
                "beta": beta,
                "required_rms_below": threshold,
                "passed": result["kappa_rms_fixed_x"] < threshold,
            }
        )
        results.append(result)
    summary = {
        "setup": {"geometry": "3x3", "U": 7, "m": 80, "R": 128, "phix": 0.25, "phiy": 0.25},
        "results": results,
        "passed": all(result["passed"] for result in results),
    }
    args.output.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    if not summary["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
