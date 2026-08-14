#!/usr/bin/python3.11
"""Print achieved-R/checkpoint/Slurm status for all sixteen twists."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import tempfile
from pathlib import Path


DEFAULT_ROOT = Path("/lustre/or-scratch24/scratch/9pm/ftlm_codex/lx4_ly4_u7_campaign")


def inspect(reducer: Path, checkpoint: Path, samples: int):
    if not checkpoint.exists() or not Path(str(checkpoint) + ".meta").exists():
        return None, "missing"
    with tempfile.NamedTemporaryFile(
        prefix="status-", suffix=".json", dir=checkpoint.parent, delete=False
    ) as handle:
        path = Path(handle.name)
    try:
        result = subprocess.run(
            [
                str(reducer), "--checkpoint", str(checkpoint), "--samples", str(samples),
                "--status-json", str(path),
            ],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        if result.returncode != 0:
            return None, "ERROR:" + (result.stderr.strip() or result.stdout.strip())
        return json.loads(path.read_text()), "ok"
    finally:
        path.unlink(missing_ok=True)


def slurm_state(name: str):
    result = subprocess.run(
        ["squeue", "-h", "-n", name, "-o", "%i:%T:%M:%R"],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    return result.stdout.strip().splitlines()[0] if result.stdout.strip() else "-"


def next_text(item):
    if not item:
        return "complete"
    suffix = "exact" if item["exact"] else f"r={item['sample_id']}"
    return f"({item['n_up']},{item['n_down']},k={item['mx']},{item['my']},{suffix})"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, required=True)
    parser.add_argument("--campaign-root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--manifest", type=Path)
    args = parser.parse_args()
    root = args.campaign_root.resolve()
    manifest = args.manifest or root / "source/campaigns/lx4_ly4_u7/twists.csv"
    reducer = root / "bin/ftlm_reduce_checkpoint"
    with manifest.open(newline="") as stream:
        rows = list(csv.DictReader(stream))

    print("twist target_R min_R progress   next_missing                         slurm")
    for row in rows:
        twist_id = row["twist_id"]
        checkpoint = root / f"runs/twist_{twist_id}/twist_{twist_id}.ftlmcp"
        status, diagnostic = inspect(reducer, checkpoint, args.samples)
        slurm = slurm_state(f"f4x4_{twist_id}")
        if status is None:
            print(f"{twist_id:>5} {args.samples:>8} {0:>5} {'0.0%':>8}   {diagnostic[:35]:<35} {slurm}")
            continue
        progress = f"{100.0 * status['weighted_progress']:.1f}%"
        print(
            f"{twist_id:>5} {args.samples:>8} {status['minimum_complete_R']:>5} "
            f"{progress:>8}   {next_text(status['next_missing']):<35} {slurm}"
        )


if __name__ == "__main__":
    main()
