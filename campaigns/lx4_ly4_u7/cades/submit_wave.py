#!/usr/bin/python3.11
"""Submit one independent CADES high-memory job for each incomplete twist."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import tempfile
from pathlib import Path


DEFAULT_ROOT = Path("/lustre/or-scratch24/scratch/9pm/ftlm_codex/lx4_ly4_u7_campaign")


def checkpoint_status(reducer: Path, checkpoint: Path, samples: int):
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
    parser.add_argument("--campaign-root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--account", default="cnms")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--allow-ungated", action="store_true")
    args = parser.parse_args()
    if args.samples <= 0:
        raise SystemExit("--samples must be positive")

    root = args.campaign_root.resolve()
    manifest = args.manifest or root / "source/campaigns/lx4_ly4_u7/twists.csv"
    reducer = root / "bin/ftlm_reduce_checkpoint"
    job_script = root / "source/campaigns/lx4_ly4_u7/cades/job_one_twist.sbatch"
    gate_path = root / "RESOURCE_GATE_PASSED.json"
    if not args.allow_ungated and not gate_path.exists():
        raise SystemExit(f"resource gate is absent: {gate_path}")
    threads = 16
    if gate_path.exists():
        gate = json.loads(gate_path.read_text())
        if not gate.get("passed"):
            raise SystemExit(f"resource gate did not pass: {gate_path}")
        threads = int(gate["production_threads"])

    root.joinpath("logs").mkdir(parents=True, exist_ok=True)
    submitted = []
    complete = []
    active = []
    with manifest.open(newline="") as stream:
        twists = list(csv.DictReader(stream))
    for row in twists:
        twist_id = row["twist_id"]
        run_dir = root / f"runs/twist_{twist_id}"
        checkpoint = run_dir / f"twist_{twist_id}.ftlmcp"
        try:
            status = checkpoint_status(reducer, checkpoint, args.samples)
        except RuntimeError as error:
            print(f"ERROR twist={twist_id} {error}")
            continue
        if status and status["complete"]:
            complete.append(twist_id)
            print(f"REUSE twist={twist_id} target_R={args.samples} state=complete")
            continue
        job_name = f"f4x4_{twist_id}"
        running = active_job(job_name)
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
                f"SEED={row['seed']}",
                f"TARGET_R={args.samples}",
                f"FTLM_THREADS={threads}",
            ]
        )
        command = [
            "sbatch", "--parsable", "--account", args.account,
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
            print(f"SUBMITTED twist={twist_id} job={job_id} target_R={args.samples}")
    print(
        f"SUMMARY target_R={args.samples} complete={len(complete)} "
        f"active={len(active)} submitted={len(submitted)} total={len(twists)}"
    )


if __name__ == "__main__":
    main()
