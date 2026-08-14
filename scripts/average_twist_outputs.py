#!/usr/bin/env python3
import argparse
import csv
import math
from collections import defaultdict


def logsumexp(values):
    top = max(values)
    return top + math.log(sum(math.exp(value - top) for value in values))


def sample_spread(values, weights):
    """Frequency-weighted spread, equivalent to expanding integer weights."""
    count = sum(weights)
    if count <= 1:
        return 0.0, 0.0
    mean = sum(weight * value for value, weight in zip(values, weights)) / count
    std = math.sqrt(
        sum(weight * (value - mean) ** 2 for value, weight in zip(values, weights))
        / (count - 1)
    )
    return std, std / math.sqrt(count)


def detect_observable(fieldnames):
    names = set(fieldnames or [])
    has_sigma = "sigma" in names
    has_sigma_dc = "sigma_dc" in names
    if has_sigma and has_sigma_dc:
        raise SystemExit("Input files cannot mix sigma and sigma_dc columns in one file")
    if has_sigma:
        return "sigma"
    if has_sigma_dc:
        return "sigma_dc"
    return None


def read_rows(paths, single_beta, weights):
    layout_seen = False
    expected_observable = None
    expected_has_omega = None
    for path, input_weight in zip(paths, weights):
        with open(path, newline="") as f:
            reader = csv.DictReader(f)
            fieldnames = reader.fieldnames or []
            has_beta = "beta" in fieldnames
            has_omega = "omega" in fieldnames
            observable = detect_observable(fieldnames)
            if not layout_seen:
                expected_observable = observable
                expected_has_omega = has_omega
                layout_seen = True
            elif observable != expected_observable or has_omega != expected_has_omega:
                raise SystemExit("All input files must have the same observable/omega column layout")
            if not has_beta and single_beta is None:
                raise SystemExit(
                    f"{path} has no beta column; pass --single-beta for legacy single-beta files"
                )
            for row in reader:
                out = {
                    "path": path,
                    "input_weight": input_weight,
                    "beta": float(row["beta"]) if has_beta else float(single_beta),
                    "mu": float(row["mu"]),
                    "omega": float(row["omega"]) if has_omega else None,
                    "n": float(row["n"]),
                    "charge_correlation": float(row["charge_correlation"]),
                    "compressibility": float(row["compressibility"]),
                    "log_partition": float(row["log_partition"])
                    if "log_partition" in row and row["log_partition"] != ""
                    else None,
                    "observable_name": observable,
                    "observable_value": float(row[observable]) if observable else None,
                }
                yield out


def main():
    parser = argparse.ArgumentParser(
        description="Average or partition-sum FTLM CSV files over twists, including n(mu), optical sigma(omega), and DC sigma."
    )
    parser.add_argument("inputs", nargs="+", help="twist-resolved CSV files")
    parser.add_argument("--output", required=True, help="output CSV path")
    parser.add_argument(
        "--mode",
        choices=["observable-average", "partition-sum"],
        default="observable-average",
        help="average observables equally, or combine twists by summing/weighting partition functions",
    )
    parser.add_argument("--sites", type=int, help="number of lattice sites; required for partition-sum")
    parser.add_argument(
        "--single-beta",
        type=float,
        help="beta value for legacy single-beta files without a beta column",
    )
    parser.add_argument(
        "--expected-twists",
        type=int,
        help="refuse output unless every group has this many independent input files",
    )
    parser.add_argument(
        "--weights",
        help="comma-separated positive integer multiplicities, one per input file",
    )
    parser.add_argument(
        "--expected-weight",
        type=int,
        help="refuse output unless the sum of input multiplicities has this value",
    )
    args = parser.parse_args()

    if args.mode == "partition-sum" and (args.sites is None or args.sites <= 0):
        raise SystemExit("--sites is required and must be positive for --mode partition-sum")

    if args.weights:
        try:
            input_weights = [int(value) for value in args.weights.split(",")]
        except ValueError as error:
            raise SystemExit("--weights must be comma-separated positive integers") from error
        if len(input_weights) != len(args.inputs):
            raise SystemExit("--weights must contain exactly one value per input file")
        if any(weight <= 0 for weight in input_weights):
            raise SystemExit("--weights values must be positive")
    else:
        input_weights = [1] * len(args.inputs)

    groups = defaultdict(list)
    observable_name = None
    has_omega = False
    for row in read_rows(args.inputs, args.single_beta, input_weights):
        observable_name = row["observable_name"] if observable_name is None else observable_name
        has_omega = row["omega"] is not None
        key = (row["beta"], row["mu"], row["omega"]) if has_omega else (row["beta"], row["mu"])
        groups[key].append(row)

    if args.expected_twists is not None and args.expected_twists <= 0:
        raise SystemExit("--expected-twists must be positive")
    if args.expected_weight is not None and args.expected_weight <= 0:
        raise SystemExit("--expected-weight must be positive")
    for key, rows in groups.items():
        paths = [row["path"] for row in rows]
        if len(set(paths)) != len(paths):
            raise SystemExit(f"duplicate input row for grid key {key}")
        if args.expected_twists is not None and len(rows) != args.expected_twists:
            raise SystemExit(
                f"grid key {key} has {len(rows)} twists; expected {args.expected_twists}"
            )
        weight_sum = sum(row["input_weight"] for row in rows)
        if args.expected_weight is not None and weight_sum != args.expected_weight:
            raise SystemExit(
                f"grid key {key} has total multiplicity {weight_sum}; "
                f"expected {args.expected_weight}"
            )

    with open(args.output, "w", newline="") as f:
        fields = ["beta", "mu"]
        if has_omega:
            fields.append("omega")
        if observable_name:
            fields.append(observable_name)
        fields += [
            "n",
            "x",
            "charge_correlation",
            "compressibility",
            "n_std",
            "n_sem",
            "x_std",
            "x_sem",
            "compressibility_std",
            "compressibility_sem",
            "twist_count",
            "twist_weight_sum",
            "mode",
            "log_partition",
        ]
        if observable_name:
            insertion = fields.index("twist_count")
            fields[insertion:insertion] = [
                f"{observable_name}_std",
                f"{observable_name}_sem",
            ]
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()

        for key, rows in sorted(groups.items()):
            beta = key[0]
            mu = key[1]
            omega = key[2] if has_omega else None
            if args.mode == "observable-average":
                weights = [row["input_weight"] for row in rows]
                weight_sum = sum(weights)
                n = sum(weight * row["n"] for weight, row in zip(weights, rows)) / weight_sum
                charge = sum(
                    weight * row["charge_correlation"]
                    for weight, row in zip(weights, rows)
                ) / weight_sum
                kappa = sum(
                    weight * row["compressibility"]
                    for weight, row in zip(weights, rows)
                ) / weight_sum
                obs = (
                    sum(
                        weight * row["observable_value"]
                        for weight, row in zip(weights, rows)
                    ) / weight_sum
                    if observable_name
                    else None
                )
                n_std, n_sem = sample_spread([row["n"] for row in rows], weights)
                kappa_std, kappa_sem = sample_spread(
                    [row["compressibility"] for row in rows], weights
                )
                if observable_name:
                    obs_std, obs_sem = sample_spread(
                        [row["observable_value"] for row in rows], weights
                    )
                weighted_logs = [
                    row["log_partition"] + math.log(row["input_weight"])
                    for row in rows
                    if row["log_partition"] is not None
                ]
                log_weight = sum(
                    row["input_weight"]
                    for row in rows
                    if row["log_partition"] is not None
                )
                log_partition = (
                    logsumexp(weighted_logs) - math.log(log_weight)
                    if weighted_logs
                    else ""
                )
            else:
                if any(row["log_partition"] is None for row in rows):
                    raise SystemExit("partition-sum mode requires log_partition in every input row")
                logz = [
                    row["log_partition"] + math.log(row["input_weight"])
                    for row in rows
                ]
                total_logz = logsumexp(logz)
                weights = [math.exp(value - total_logz) for value in logz]
                mean_n_total = 0.0
                mean_n2_total = 0.0
                obs = 0.0 if observable_name else None
                for weight, row in zip(weights, rows):
                    particles = row["n"] * args.sites
                    particles2 = row["charge_correlation"] * args.sites + particles * particles
                    mean_n_total += weight * particles
                    mean_n2_total += weight * particles2
                    if observable_name:
                        obs += weight * row["observable_value"]
                n = mean_n_total / args.sites
                charge = max(0.0, (mean_n2_total - mean_n_total * mean_n_total) / args.sites)
                kappa = beta * charge
                log_partition = total_logz
                n_std = n_sem = kappa_std = kappa_sem = ""
                if observable_name:
                    obs_std = obs_sem = ""

            out = {
                "beta": beta,
                "mu": mu,
                "n": n,
                "x": 1.0 - n,
                "charge_correlation": charge,
                "compressibility": kappa,
                "n_std": n_std,
                "n_sem": n_sem,
                "x_std": n_std,
                "x_sem": n_sem,
                "compressibility_std": kappa_std,
                "compressibility_sem": kappa_sem,
                "twist_count": len(rows),
                "twist_weight_sum": sum(row["input_weight"] for row in rows),
                "mode": args.mode,
                "log_partition": log_partition,
            }
            if has_omega:
                out["omega"] = omega
            if observable_name:
                out[observable_name] = obs
                out[f"{observable_name}_std"] = obs_std
                out[f"{observable_name}_sem"] = obs_sem
            writer.writerow(out)


if __name__ == "__main__":
    main()
