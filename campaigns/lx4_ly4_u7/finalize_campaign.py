#!/usr/bin/env python3
"""Validate, symmetry-weight, and plot a complete target-R campaign."""

from __future__ import annotations

import argparse
import csv
import json
import math
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


def validate_manifest(twists, quadrature):
    if not twists:
        raise RuntimeError("twist manifest is empty")
    lx = int(quadrature["cluster_lx"])
    ly = int(quadrature["cluster_ly"])
    if lx != 4 or ly != 4:
        raise RuntimeError(f"expected a 4x4 quadrature, found {lx}x{ly}")
    grid_size = int(quadrature["grid_nx"])
    if grid_size != 20 or int(quadrature["grid_ny"]) != grid_size:
        raise RuntimeError("expected a 20x20 reduced-BZ grid")
    if len(twists) != int(quadrature["representatives"]):
        raise RuntimeError("manifest representative count disagrees with metadata")
    ids = set()
    total_weight = 0
    for row in twists:
        twist_id = row["twist_id"]
        if twist_id in ids:
            raise RuntimeError(f"duplicate twist ID: {twist_id}")
        ids.add(twist_id)
        ix = int(row["grid_ix"])
        iy = int(row["grid_iy"])
        phix = float(row["phix"])
        phiy = float(row["phiy"])
        kx_over_pi = float(row["kx_over_pi"])
        ky_over_pi = float(row["ky_over_pi"])
        if not (0 <= ix <= iy <= grid_size // 2):
            raise RuntimeError(
                f"twist {twist_id} lies outside 0 <= kx <= ky <= pi/4: {row}"
            )
        if not math.isclose(phix, ix / grid_size, abs_tol=1e-12):
            raise RuntimeError(f"twist {twist_id} has inconsistent phix/grid_ix")
        if not math.isclose(phiy, iy / grid_size, abs_tol=1e-12):
            raise RuntimeError(f"twist {twist_id} has inconsistent phiy/grid_iy")
        if not math.isclose(kx_over_pi, 2.0 * phix / lx, abs_tol=1e-12):
            raise RuntimeError(f"twist {twist_id} has inconsistent kx/phix")
        if not math.isclose(ky_over_pi, 2.0 * phiy / ly, abs_tol=1e-12):
            raise RuntimeError(f"twist {twist_id} has inconsistent ky/phiy")
        orbit = set()
        for sx in (-1, 1):
            for sy in (-1, 1):
                orbit.add(((sx * ix) % grid_size, (sy * iy) % grid_size))
                orbit.add(((sx * iy) % grid_size, (sy * ix) % grid_size))
        expected_multiplicity = len(orbit)
        multiplicity = int(row["multiplicity"])
        if multiplicity != expected_multiplicity:
            raise RuntimeError(
                f"twist {twist_id} has multiplicity {multiplicity}, "
                f"expected {expected_multiplicity}"
            )
        total_weight += multiplicity
    expected_weight = int(quadrature["effective_twists"])
    if total_weight != expected_weight:
        raise RuntimeError(
            f"twist multiplicities sum to {total_weight}, expected {expected_weight}"
        )
    return total_weight


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--campaign-root", type=Path, required=True)
    parser.add_argument("--samples", type=int, required=True)
    parser.add_argument("--lanczos-steps", type=int, default=80)
    parser.add_argument("--checkpoint-series", default="legacy")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--quadrature-metadata", type=Path)
    parser.add_argument("--skip-plot", action="store_true")
    args = parser.parse_args()
    root = args.campaign_root.resolve()
    tag = f"{args.samples:03d}"
    m_tag = f"{args.lanczos_steps:03d}"
    result_tag = (
        f"R{tag}"
        if args.checkpoint_series == "legacy"
        else f"m{m_tag}_R{tag}"
    )
    campaign_dir = root / "source/campaigns/lx4_ly4_u7"
    manifest = (args.manifest or campaign_dir / "twists.csv").resolve()
    metadata_path = (
        args.quadrature_metadata or campaign_dir / "twist_quadrature.json"
    ).resolve()
    with manifest.open(newline="") as stream:
        twists = list(csv.DictReader(stream))
    quadrature = json.loads(metadata_path.read_text())
    try:
        effective_twists = validate_manifest(twists, quadrature)
    except RuntimeError as error:
        raise SystemExit(f"invalid reduced-BZ quadrature: {error}") from error
    inputs = [
        root / (
            f"runs/twist_{row['twist_id']}/twist_{row['twist_id']}_thermo_"
            f"{result_tag}.csv"
        )
        for row in twists
    ]
    weights = [int(row.get("multiplicity", "1")) for row in twists]
    missing = [str(path) for path in inputs if not path.exists()]
    if missing:
        raise SystemExit("missing twist outputs:\n" + "\n".join(missing))
    for path in inputs:
        validate_twist(path)

    source = root / "source"
    output_dir = (args.output_dir or root / "averages").resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    quadrature_tag = quadrature["tag"]
    average = output_dir / f"twist_average_{quadrature_tag}_{result_tag}.csv"
    subprocess.run(
        [
            sys.executable,
            str(source / "scripts/average_twist_outputs.py"),
            *map(str, inputs),
            "--mode", "observable-average",
            "--weights", ",".join(map(str, weights)),
            "--expected-twists", str(len(inputs)),
            "--expected-weight", str(effective_twists),
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
                "--output",
                str(output_dir / f"compressibility_vs_x_{quadrature_tag}_{result_tag}.png"),
            ],
            check=True,
        )
    print(
        f"quadrature={quadrature_tag} validated_representatives={len(inputs)} "
        f"effective_twists={effective_twists} average={average}"
    )


if __name__ == "__main__":
    main()
