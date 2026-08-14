#!/usr/bin/env python3
"""Plot twist-averaged compressibility against x=1-n."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt


def read_groups(path: Path, expected_mu_count: int | None):
    groups = defaultdict(list)
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        required = {"beta", "n", "compressibility"}
        if not required.issubset(reader.fieldnames or []):
            raise SystemExit(f"{path} lacks required columns: {sorted(required)}")
        for row in reader:
            beta = float(row["beta"])
            groups[beta].append(
                {
                    "x": float(row["x"]) if row.get("x") else 1.0 - float(row["n"]),
                    "kappa": float(row["compressibility"]),
                    "std": float(row["compressibility_std"])
                    if row.get("compressibility_std") else None,
                    "sem": float(row["compressibility_sem"])
                    if row.get("compressibility_sem") else None,
                }
            )
    if expected_mu_count is not None:
        for beta, rows in groups.items():
            if len(rows) != expected_mu_count:
                raise SystemExit(
                    f"{path}: beta={beta:g} has {len(rows)} rows; expected {expected_mu_count}"
                )
    return groups


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="+", type=Path, help="twist-average CSV files")
    parser.add_argument("--labels", nargs="*", help="optional dataset labels, e.g. R=8 R=16")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--x-min", type=float, default=0.0)
    parser.add_argument("--x-max", type=float, default=0.35)
    parser.add_argument("--expected-mu-count", type=int, default=281)
    parser.add_argument("--error", choices=["none", "std", "sem"], default="sem")
    args = parser.parse_args()
    if args.labels and len(args.labels) != len(args.inputs):
        raise SystemExit("--labels must contain one label per input")

    labels = args.labels or ([""] if len(args.inputs) == 1 else [path.stem for path in args.inputs])
    all_groups = [read_groups(path, args.expected_mu_count) for path in args.inputs]
    beta_values = sorted({beta for groups in all_groups for beta in groups})
    colors = plt.get_cmap("tab10")
    color_by_beta = {beta: colors(index) for index, beta in enumerate(beta_values)}
    line_styles = ["-", "--", ":", "-."]

    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.linewidth": 0.9,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )
    fig, ax = plt.subplots(figsize=(3.45, 2.8), constrained_layout=True)
    for dataset_index, (groups, dataset_label) in enumerate(zip(all_groups, labels)):
        for beta in sorted(groups):
            rows = sorted(
                (row for row in groups[beta] if args.x_min <= row["x"] <= args.x_max),
                key=lambda row: row["x"],
            )
            if not rows:
                continue
            x = [row["x"] for row in rows]
            y = [row["kappa"] for row in rows]
            temperature = 1.0 / beta
            legend = rf"$T/t={temperature:.2f}$"
            if dataset_label:
                legend += f" ({dataset_label})"
            ax.plot(
                x,
                y,
                color=color_by_beta[beta],
                linestyle=line_styles[dataset_index % len(line_styles)],
                linewidth=1.45,
                label=legend,
            )
            if args.error != "none" and all(row[args.error] is not None for row in rows):
                error = [row[args.error] for row in rows]
                lower = [value - spread for value, spread in zip(y, error)]
                upper = [value + spread for value, spread in zip(y, error)]
                ax.fill_between(x, lower, upper, color=color_by_beta[beta], alpha=0.14, linewidth=0)

    ax.set_xlim(args.x_min, args.x_max)
    ax.set_xlabel(r"Hole doping $x=1-n$")
    ax.set_ylabel(r"Compressibility $\kappa$")
    ax.legend(frameon=False, fontsize=8.5)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=300)
    if args.output.suffix.lower() == ".png":
        fig.savefig(args.output.with_suffix(".pdf"))
    print(args.output.resolve())


if __name__ == "__main__":
    main()
