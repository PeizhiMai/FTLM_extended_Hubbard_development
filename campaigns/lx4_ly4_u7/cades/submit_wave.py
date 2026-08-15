#!/usr/bin/python3.11
"""Submit one independent CADES high-memory job for each incomplete twist."""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import subprocess
import tempfile
from pathlib import Path


DEFAULT_ROOT = Path("/lustre/or-scratch24/scratch/9pm/ftlm_codex/lx4_ly4_u7_campaign")
CADES_HIGH_MEM_CORES = 36


def validate_rbz_manifest(rows, metadata):
    if metadata.get("tag") != "rbz20_gamma_dk_pi40":
        raise RuntimeError(f"unexpected quadrature tag: {metadata.get('tag')}")
    if len(rows) != int(metadata["representatives"]):
        raise RuntimeError("manifest representative count disagrees with metadata")
    grid_size = int(metadata["grid_nx"])
    if grid_size != 20 or int(metadata["grid_ny"]) != grid_size:
        raise RuntimeError("expected a 20x20 reduced-BZ grid")
    ids = set()
    weight = 0
    for row in rows:
        twist_id = row["twist_id"]
        if twist_id in ids:
            raise RuntimeError(f"duplicate twist ID: {twist_id}")
        ids.add(twist_id)
        ix = int(row["grid_ix"])
        iy = int(row["grid_iy"])
        phix = float(row["phix"])
        phiy = float(row["phiy"])
        if not (0 <= ix <= iy <= grid_size // 2):
            raise RuntimeError(
                f"twist {twist_id} is outside 0 <= kx <= ky <= pi/4"
            )
        if not math.isclose(phix, ix / grid_size, abs_tol=1e-12):
            raise RuntimeError(f"twist {twist_id} has inconsistent phix/grid_ix")
        if not math.isclose(phiy, iy / grid_size, abs_tol=1e-12):
            raise RuntimeError(f"twist {twist_id} has inconsistent phiy/grid_iy")
        if not math.isclose(float(row["kx_over_pi"]), ix / 40.0, abs_tol=1e-12):
            raise RuntimeError(f"twist {twist_id} has inconsistent kx/grid_ix")
        if not math.isclose(float(row["ky_over_pi"]), iy / 40.0, abs_tol=1e-12):
            raise RuntimeError(f"twist {twist_id} has inconsistent ky/grid_iy")
        orbit = set()
        for sx in (-1, 1):
            for sy in (-1, 1):
                orbit.add(((sx * ix) % grid_size, (sy * iy) % grid_size))
                orbit.add(((sx * iy) % grid_size, (sy * ix) % grid_size))
        expected = len(orbit)
        if int(row["multiplicity"]) != expected:
            raise RuntimeError(
                f"twist {twist_id} has multiplicity {row['multiplicity']}, "
                f"expected {expected}"
            )
        weight += expected
    if weight != int(metadata["effective_twists"]):
        raise RuntimeError(
            f"manifest weight {weight} != {metadata['effective_twists']}"
        )
    return weight


def checkpoint_status(reducer: Path, checkpoint: Path, samples: int, lanczos_steps: int):
    if not checkpoint.exists() or not Path(str(checkpoint) + ".meta").exists():
        return None
    checkpoint.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        prefix="status-", suffix=".json", dir=checkpoint.parent, delete=False
    ) as handle:
        status_path = Path(handle.name)
    try:
        result = subprocess.run(
            [
                str(reducer), "--checkpoint", str(checkpoint), "--samples", str(samples),
                "--lanczos-steps", str(lanczos_steps),
                "--status-json", str(status_path),
            ],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"cannot inspect {checkpoint}: {result.stderr.strip() or result.stdout.strip()}"
            )
        return json.loads(status_path.read_text())
    finally:
        status_path.unlink(missing_ok=True)


def active_job(job_name: str) -> str | None:
    result = subprocess.run(
        ["squeue", "-h", "-n", job_name, "-o", "%i|%T|%R"],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    line = result.stdout.strip().splitlines()
    return line[0] if line else None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, required=True, help="target total R")
    parser.add_argument(
        "--threads",
        type=int,
        help=(
            "OpenMP worker count; defaults to the resource-gate recommendation. "
            "The Slurm CPU request is capped at the 36 physical high-memory cores."
        ),
    )
    parser.add_argument("--lanczos-max-steps", type=int, default=80)
    parser.add_argument(
        "--lanczos-save-steps",
        help="comma-separated prefixes saved from one recurrence (for example 80,120,300)",
    )
    parser.add_argument(
        "--checkpoint-series",
        default="legacy",
        help="separate checkpoint namespace; use mext for multi-prefix runs",
    )
    parser.add_argument(
        "--only-block",
        help="optional resource-probe block NUP,NDOWN,MX,MY",
    )
    parser.add_argument("--campaign-root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--account", default="cnms")
    parser.add_argument(
        "--twist-id",
        action="append",
        help="submit only this three-digit twist ID (repeatable); default is all twists",
    )
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--allow-ungated", action="store_true")
    args = parser.parse_args()
    if args.samples <= 0:
        raise SystemExit("--samples must be positive")
    if args.threads is not None and args.threads <= 0:
        raise SystemExit("--threads must be positive")
    if args.lanczos_max_steps <= 0:
        raise SystemExit("--lanczos-max-steps must be positive")
    if not args.checkpoint_series.replace("_", "").replace("-", "").isalnum():
        raise SystemExit("--checkpoint-series may contain only letters, digits, '_' and '-'")
    save_steps = args.lanczos_save_steps
    if args.checkpoint_series != "legacy" and not save_steps:
        save_steps = str(args.lanczos_max_steps)
    if save_steps:
        try:
            parsed_steps = sorted({int(value) for value in save_steps.split(",")})
        except ValueError as error:
            raise SystemExit("--lanczos-save-steps must be comma-separated integers") from error
        if not parsed_steps or parsed_steps[0] <= 0 or parsed_steps[-1] > args.lanczos_max_steps:
            raise SystemExit("saved Lanczos steps must be in [1,--lanczos-max-steps]")
        if args.lanczos_max_steps not in parsed_steps:
            parsed_steps.append(args.lanczos_max_steps)
            parsed_steps.sort()
        save_steps = ",".join(map(str, parsed_steps))
    only_block = None
    if args.only_block:
        try:
            only_block = tuple(int(value) for value in args.only_block.split(","))
        except ValueError as error:
            raise SystemExit("--only-block expects NUP,NDOWN,MX,MY") from error
        if len(only_block) != 4:
            raise SystemExit("--only-block expects NUP,NDOWN,MX,MY")

    root = args.campaign_root.resolve()
    manifest = args.manifest or root / "source/campaigns/lx4_ly4_u7/twists.csv"
    reducer = root / "bin/ftlm_reduce_checkpoint"
    job_script = root / "source/campaigns/lx4_ly4_u7/cades/job_one_twist.sbatch"
    gate_path = root / (
        "RESOURCE_GATE_PASSED.json"
        if args.checkpoint_series == "legacy"
        else "RESOURCE_GATE_M300_PASSED.json"
    )
    if not args.allow_ungated and not gate_path.exists():
        raise SystemExit(f"resource gate is absent: {gate_path}")
    threads = 16
    if gate_path.exists():
        gate = json.loads(gate_path.read_text())
        if not gate.get("passed"):
            raise SystemExit(f"resource gate did not pass: {gate_path}")
        threads = int(gate["production_threads"])
    if args.threads is not None:
        threads = args.threads
    allocated_cpus = min(threads, CADES_HIGH_MEM_CORES)

    root.joinpath("logs").mkdir(parents=True, exist_ok=True)
    submitted = []
    complete = []
    active = []
    with manifest.open(newline="") as stream:
        twists = list(csv.DictReader(stream))
    metadata_path = manifest.with_name("twist_quadrature.json")
    if not metadata_path.exists():
        raise SystemExit(f"reduced-BZ quadrature metadata is absent: {metadata_path}")
    metadata = json.loads(metadata_path.read_text())
    try:
        effective_twists = validate_rbz_manifest(twists, metadata)
    except (KeyError, TypeError, ValueError, RuntimeError) as error:
        raise SystemExit(f"invalid reduced-BZ manifest: {error}") from error
    if args.twist_id:
        requested = {value.zfill(3) for value in args.twist_id}
        known = {row["twist_id"] for row in twists}
        unknown = requested - known
        if unknown:
            raise SystemExit(f"unknown twist IDs: {sorted(unknown)}")
        twists = [row for row in twists if row["twist_id"] in requested]
    selected_weight = sum(int(row["multiplicity"]) for row in twists)
    for row in twists:
        twist_id = row["twist_id"]
        run_dir = root / f"runs/twist_{twist_id}"
        checkpoint_suffix = "" if args.checkpoint_series == "legacy" else f"_{args.checkpoint_series}"
        checkpoint = run_dir / f"twist_{twist_id}{checkpoint_suffix}.ftlmcp"
        try:
            status = checkpoint_status(
                reducer, checkpoint, args.samples, args.lanczos_max_steps
            )
        except RuntimeError as error:
            print(f"ERROR twist={twist_id} {error}")
            continue
        selected_complete = False
        if status and only_block:
            selected_complete = all(
                (
                    item["n_up"], item["n_down"], item["mx"], item["my"]
                ) != only_block
                for item in status["missing_blocks"]
            )
        if status and (status["complete"] or selected_complete):
            complete.append(twist_id)
            scope = "selected-block" if selected_complete and not status["complete"] else "run"
            print(
                f"REUSE twist={twist_id} target_R={args.samples} "
                f"target_m={args.lanczos_max_steps} scope={scope} state=complete"
            )
            continue
        job_name = (
            f"f4x4_{twist_id}"
            if args.checkpoint_series == "legacy"
            else f"f4x4_{args.checkpoint_series}_{twist_id}"
        )
        running = None if args.dry_run and shutil.which("squeue") is None else active_job(job_name)
        if running:
            active.append((twist_id, running))
            print(f"ACTIVE twist={twist_id} job={running}")
            continue
        exports = ",".join(
            [
                "ALL",
                f"CAMPAIGN_ROOT={root}",
                f"TWIST_ID={twist_id}",
                f"PHIX={row['phix']}",
                f"PHIY={row['phiy']}",
                f"TWIST_KX_OVER_PI={row['kx_over_pi']}",
                f"TWIST_KY_OVER_PI={row['ky_over_pi']}",
                f"TWIST_MULTIPLICITY={row['multiplicity']}",
                f"TWIST_QUADRATURE_TAG={metadata['tag']}",
                f"SEED={row['seed']}",
                f"TARGET_R={args.samples}",
                f"FTLM_THREADS={threads}",
                f"LANCZOS_MAX_STEPS={args.lanczos_max_steps}",
                f"LANCZOS_SAVE_STEPS={(save_steps or '').replace(',', ':')}",
                f"CHECKPOINT_SERIES={args.checkpoint_series}",
            ]
        )
        if only_block:
            exports += "," + ",".join(
                [
                    f"ONLY_BLOCK_NUP={only_block[0]}",
                    f"ONLY_BLOCK_NDOWN={only_block[1]}",
                    f"ONLY_BLOCK_MX={only_block[2]}",
                    f"ONLY_BLOCK_MY={only_block[3]}",
                ]
            )
        command = [
            "sbatch", "--parsable", "--account", args.account,
            "--cpus-per-task", str(allocated_cpus),
            "--job-name", job_name,
            "--output", str(root / "logs/%x-%j.out"),
            "--error", str(root / "logs/%x-%j.err"),
            "--export", exports,
            str(job_script),
        ]
        if args.dry_run:
            print("DRY_RUN", " ".join(command))
            submitted.append((twist_id, "dry-run"))
        else:
            result = subprocess.run(command, text=True, stdout=subprocess.PIPE, check=True)
            job_id = result.stdout.strip().split(";")[0]
            submitted.append((twist_id, job_id))
            print(
                f"SUBMITTED twist={twist_id} job={job_id} target_R={args.samples} "
                f"target_m={args.lanczos_max_steps} series={args.checkpoint_series} "
                f"threads={threads} allocated_cpus={allocated_cpus} mem=250G"
            )
    print(
        f"SUMMARY target_R={args.samples} target_m={args.lanczos_max_steps} "
        f"series={args.checkpoint_series} threads={threads} "
        f"allocated_cpus={allocated_cpus} complete={len(complete)} "
        f"active={len(active)} submitted={len(submitted)} total={len(twists)} "
        f"quadrature={metadata['tag']} selected_weight={selected_weight} "
        f"full_effective_twists={effective_twists}"
    )


if __name__ == "__main__":
    main()
