#!/usr/bin/env python3
import argparse
import csv
import math
from collections import defaultdict


def logsumexp(values):
    top = max(values)
    return top + math.log(sum(math.exp(value - top) for value in values))


def read_rows(paths, single_beta):
    for path in paths:
        with open(path, newline="") as f:
            reader = csv.DictReader(f)
            has_beta = "beta" in (reader.fieldnames or [])
            if not has_beta and single_beta is None:
                raise SystemExit(
                    f"{path} has no beta column; pass --single-beta for legacy single-beta files"
                )
            for row in reader:
                beta = float(row["beta"]) if has_beta else float(single_beta)
                yield {
                    "path": path,
                    "beta": beta,
                    "mu": float(row["mu"]),
                    "n": float(row["n"]),
                    "charge_correlation": float(row["charge_correlation"]),
                    "compressibility": float(row["compressibility"]),
                    "log_partition": float(row["log_partition"])
                    if "log_partition" in row and row["log_partition"] != ""
                    else None,
                }


def main():
    parser = argparse.ArgumentParser(
        description="Average or partition-sum FTLM n(mu) CSV files over twists."
    )
    parser.add_argument("inputs", nargs="+", help="twist-resolved CSV files")
    parser.add_argument("--output", required=True, help="output CSV path")
    parser.add_argument(
        "--mode",
        choices=["observable-average", "partition-sum"],
        default="observable-average",
        help="average observables equally, or combine twists by summing partition functions",
    )
    parser.add_argument("--sites", type=int, help="number of lattice sites; required for partition-sum")
    parser.add_argument(
        "--single-beta",
        type=float,
        help="beta value for legacy single-beta files without a beta column",
    )
    args = parser.parse_args()

    if args.mode == "partition-sum" and (args.sites is None or args.sites <= 0):
        raise SystemExit("--sites is required and must be positive for --mode partition-sum")

    groups = defaultdict(list)
    for row in read_rows(args.inputs, args.single_beta):
        groups[(row["beta"], row["mu"])].append(row)

    with open(args.output, "w", newline="") as f:
        fields = [
            "beta",
            "mu",
            "n",
            "charge_correlation",
            "compressibility",
            "twist_count",
            "mode",
            "log_partition",
        ]
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()

        for (beta, mu), rows in sorted(groups.items()):
            if args.mode == "observable-average":
                count = len(rows)
                n = sum(row["n"] for row in rows) / count
                charge = sum(row["charge_correlation"] for row in rows) / count
                kappa = sum(row["compressibility"] for row in rows) / count
                logs = [row["log_partition"] for row in rows if row["log_partition"] is not None]
                log_partition = logsumexp(logs) - math.log(len(logs)) if logs else ""
            else:
                if any(row["log_partition"] is None for row in rows):
                    raise SystemExit("partition-sum mode requires log_partition in every input row")
                logz = [row["log_partition"] for row in rows]
                total_logz = logsumexp(logz)
                weights = [math.exp(value - total_logz) for value in logz]
                mean_n_total = 0.0
                mean_n2_total = 0.0
                for weight, row in zip(weights, rows):
                    particles = row["n"] * args.sites
                    particles2 = row["charge_correlation"] * args.sites + particles * particles
                    mean_n_total += weight * particles
                    mean_n2_total += weight * particles2
                n = mean_n_total / args.sites
                charge = max(0.0, (mean_n2_total - mean_n_total * mean_n_total) / args.sites)
                kappa = beta * charge
                log_partition = total_logz

            writer.writerow(
                {
                    "beta": beta,
                    "mu": mu,
                    "n": n,
                    "charge_correlation": charge,
                    "compressibility": kappa,
                    "twist_count": len(rows),
                    "mode": args.mode,
                    "log_partition": log_partition,
                }
            )


if __name__ == "__main__":
    main()
