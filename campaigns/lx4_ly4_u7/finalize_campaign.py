#!/usr/bin/env python3
"""Validate, symmetry-weight, and plot a complete target-R campaign."""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


def validate_twist(path: Path):
    groups = defaultdict(list)
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            groups[float(row["beta"])].append(row)
    if len(groups) != 2 or any(len(rows) != 281 for rows in groups.values()):
        raise RuntimeError(f"{path}: expected two betas and 281 mu rows per beta")
    for beta, rows in groups.items():
        rows.sort(key=lambda row: float(row["mu"]))
        densities = [float(row["n"]) for row in rows]
        kappas = [float(row["compressibility"]) for row in rows]
        if any(right + 1e-9 < left for left, right in zip(densities, densities[1:])):
            raise RuntimeError(f"{path}: density is not monotone at beta={beta}")
        if min(kappas) < -1e-10:
            raise RuntimeError(f"{path}: negative compressibility at beta={beta}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--campaign-root", type=Path, required=True)
    parser.add_argument("--samples", type=int, required=True)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--skip-plot", action="store_true")
    args = parser.parse_args()
    root = args.campaign_root.resolve()
    tag = f"{args.samples:03d}"
    manifest = root / "source/campaigns/lx4_ly4_u7/twists.csv"
    with manifest.open(newline="") as stream:
        twists = list(csv.DictReader(stream))
    if not twists:
        raise SystemExit(f"empty twist manifest: {manifest}")
    for row in twists:
        if float(row["phiy"]) < float(row["phix"]):
            raise SystemExit(f"manifest contains non-representative twist: {row}")
    inputs = [
        root / f"runs/twist_{row['twist_id']}/twist_{row['twist_id']}_thermo_R{tag}.csv"
        for row in twists
    ]
    weights = [int(row.get("multiplicity", "1")) for row in twists]
    if sum(weights) != 16:
        raise SystemExit(
            f"twist symmetry multiplicities sum to {sum(weights)}, expected full-grid weight 16"
        )
    missing = [str(path) for path in inputs if not path.exists()]
    if missing:
        raise SystemExit("missing twist outputs:\n" + "\n".join(missing))
    for path in inputs:
        validate_twist(path)

    source = root / "source"
    output_dir = (args.output_dir or root / "averages").resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    average = output_dir / f"twist_average_R{tag}.csv"
    subprocess.run(
        [
            sys.executable,
            str(source / "scripts/average_twist_outputs.py"),
            *map(str, inputs),
            "--mode", "observable-average",
            "--weights", ",".join(map(str, weights)),
            "--expected-twists", str(len(inputs)),
            "--expected-weight", "16",
            "--output", str(average),
        ],
        check=True,
    )
    with average.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    for beta in {float(row["beta"]) for row in rows}:
        x = [float(row["x"]) for row in rows if float(row["beta"]) == beta]
        if min(x) > 0.0 + 2e-3 or max(x) < 0.35:
            raise RuntimeError(f"twist average lacks x=[0,0.35] coverage at beta={beta}")
    if not args.skip_plot:
        subprocess.run(
            [
                sys.executable,
                str(source / "scripts/plot_compressibility_vs_hole_doping.py"),
                str(average),
                "--output", str(output_dir / f"compressibility_vs_x_R{tag}.png"),
            ],
            check=True,
        )
    print(
        f"validated_representatives={len(inputs)} effective_twists={sum(weights)} "
        f"average={average}"
    )


if __name__ == "__main__":
    main()
