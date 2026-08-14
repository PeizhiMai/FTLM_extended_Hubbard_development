#!/usr/bin/python3.11
"""Validate one completed twist-resolved thermodynamic CSV."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path


EXPECTED_BETAS = (2.857142857142857, 12.5)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    groups = defaultdict(list)
    with args.input.open(newline="") as stream:
        reader = csv.DictReader(stream)
        required = {"beta", "mu", "n", "charge_correlation", "compressibility", "log_partition"}
        if not required.issubset(reader.fieldnames or []):
            raise SystemExit(f"missing required columns: {sorted(required - set(reader.fieldnames or []))}")
        for row in reader:
            groups[float(row["beta"])].append({key: float(value) for key, value in row.items()})

    results = []
    reasons = []
    if len(groups) != 2:
        reasons.append(f"found {len(groups)} betas instead of 2")
    for expected_beta in EXPECTED_BETAS:
        match = next((beta for beta in groups if abs(beta - expected_beta) < 1e-10), None)
        if match is None:
            reasons.append(f"missing beta={expected_beta}")
            continue
        rows = sorted(groups[match], key=lambda row: row["mu"])
        densities = [row["n"] for row in rows]
        kappas = [row["compressibility"] for row in rows]
        x = [1.0 - density for density in densities]
        finite = all(math.isfinite(value) for row in rows for value in row.values())
        monotone = all(right + 1e-9 >= left for left, right in zip(densities, densities[1:]))
        nonnegative = min(kappas, default=-1.0) >= -1e-10
        identity_error = max(
            (abs(row["compressibility"] - match * row["charge_correlation"]) for row in rows),
            default=math.inf,
        )
        coverage = min(x, default=math.inf) <= 0.002 and max(x, default=-math.inf) >= 0.35
        result = {
            "beta": match,
            "temperature": 1.0 / match,
            "rows": len(rows),
            "mu_min": rows[0]["mu"] if rows else None,
            "mu_max": rows[-1]["mu"] if rows else None,
            "n_min": min(densities, default=None),
            "n_max": max(densities, default=None),
            "kappa_min": min(kappas, default=None),
            "kappa_max": max(kappas, default=None),
            "x_min": min(x, default=None),
            "x_max": max(x, default=None),
            "finite": finite,
            "density_monotone": monotone,
            "compressibility_nonnegative": nonnegative,
            "kappa_beta_charge_max_error": identity_error,
            "covers_x_0_to_0p35": coverage,
        }
        results.append(result)
        if len(rows) != 281:
            reasons.append(f"beta={match} has {len(rows)} rows instead of 281")
        if not finite:
            reasons.append(f"beta={match} contains nonfinite values")
        if not monotone:
            reasons.append(f"beta={match} density is not monotone")
        if not nonnegative:
            reasons.append(f"beta={match} has negative compressibility")
        if identity_error > 1e-10:
            reasons.append(f"beta={match} violates kappa=beta*charge")
        if not coverage:
            reasons.append(f"beta={match} does not cover x=[0,0.35]")

    summary = {"input": str(args.input.resolve()), "passed": not reasons, "reasons": reasons, "curves": results}
    text = json.dumps(summary, indent=2) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")
    if reasons:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
