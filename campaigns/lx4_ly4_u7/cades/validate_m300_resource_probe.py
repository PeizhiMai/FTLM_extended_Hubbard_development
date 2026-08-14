#!/usr/bin/python3.11
"""Validate the separate m_max=300 largest-block resource probe."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import tempfile
from pathlib import Path


DEFAULT_ROOT = Path("/lustre/or-scratch24/scratch/9pm/ftlm_codex/lx4_ly4_u7_campaign")
EXPECTED_DIMENSION = 10_353_252
EXPECTED_SAMPLES = set(range(16))
SAVE_STEPS = (80, 120, 160, 200, 250, 300)
TARGET_BLOCK = (8, 8, 0, 0)


def atomic_json(path: Path, value):
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def block_complete_at_depth(reducer: Path, checkpoint: Path, steps: int) -> bool:
    with tempfile.NamedTemporaryFile(
        prefix=f"m{steps}-status-", suffix=".json", dir=checkpoint.parent, delete=False
    ) as handle:
        status_path = Path(handle.name)
    try:
        result = subprocess.run(
            [
                str(reducer), "--checkpoint", str(checkpoint), "--samples", "16",
                "--lanczos-steps", str(steps), "--status-json", str(status_path),
            ],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        if result.returncode != 0:
            raise RuntimeError(result.stderr.strip() or result.stdout.strip())
        status = json.loads(status_path.read_text())
        for item in status["missing_blocks"]:
            key = (item["n_up"], item["n_down"], item["mx"], item["my"])
            if key == TARGET_BLOCK:
                return False
        return True
    finally:
        status_path.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--campaign-root", type=Path, default=DEFAULT_ROOT)
    args = parser.parse_args()
    root = args.campaign_root.resolve()
    run_dir = root / "runs/twist_005"
    progress = run_dir / "progress_mext.jsonl"
    checkpoint = run_dir / "twist_005_mext.ftlmcp"
    reducer = root / "bin/ftlm_reduce_checkpoint"
    if not progress.exists():
        raise SystemExit(f"missing progress log: {progress}")
    if not checkpoint.exists() or not Path(str(checkpoint) + ".meta").exists():
        raise SystemExit(f"missing mext checkpoint or metadata: {checkpoint}")

    dimensions = set()
    peak_rss = 0.0
    sample_durations = {}
    malformed = 0
    for line in progress.read_text().splitlines():
        try:
            event = json.loads(line)
        except json.JSONDecodeError:
            malformed += 1
            continue
        if int(event.get("target_m", 0)) != 300:
            continue
        peak_rss = max(peak_rss, float(event.get("max_rss_gb", 0.0)))
        block = event.get("current_block")
        if block and (
            block.get("n_up"), block.get("n_down"), block.get("mx"), block.get("my")
        ) == TARGET_BLOCK:
            dimensions.add(int(block.get("basis_dim", -1)))
            if event.get("event") == "CHECKPOINTED" and event.get("checkpointed_sample") is not None:
                duration = float(event.get("last_sample_duration_seconds", -1.0))
                if duration >= 0.0:
                    sample_durations[int(event["checkpointed_sample"])] = duration

    reasons = []
    if dimensions != {EXPECTED_DIMENSION}:
        reasons.append(f"largest-block dimension observations were {sorted(dimensions)}")
    missing_samples = sorted(EXPECTED_SAMPLES - set(sample_durations))
    if missing_samples:
        reasons.append(f"missing durable m=300 timing for sample IDs {missing_samples}")
    if peak_rss >= 245.0:
        reasons.append(f"peak RSS {peak_rss:.2f} GiB leaves less than 5 GiB headroom")

    depth_complete = {}
    for steps in SAVE_STEPS:
        try:
            depth_complete[steps] = block_complete_at_depth(reducer, checkpoint, steps)
        except RuntimeError as error:
            depth_complete[steps] = False
            reasons.append(f"m={steps} checkpoint status failed: {error}")
    incomplete_depths = [steps for steps, complete in depth_complete.items() if not complete]
    if incomplete_depths:
        reasons.append(f"largest block lacks R=16 records at m={incomplete_depths}")

    report = {
        "passed": not reasons,
        "checkpoint_series": "mext",
        "requested_memory_gb": 250,
        "lanczos_max_steps": 300,
        "lanczos_save_steps": list(SAVE_STEPS),
        "expected_dimension": EXPECTED_DIMENSION,
        "observed_dimensions": sorted(dimensions),
        "peak_rss_gb": peak_rss,
        "sample_durations_seconds": sample_durations,
        "max_sample_seconds": max(sample_durations.values(), default=None),
        "missing_sample_ids": missing_samples,
        "depth_complete": depth_complete,
        "production_threads": 16,
        "malformed_progress_lines": malformed,
        "reasons": reasons,
    }
    report_path = root / "resource_probe_m300_report.json"
    atomic_json(report_path, report)
    gate_path = root / "RESOURCE_GATE_M300_PASSED.json"
    if report["passed"]:
        atomic_json(gate_path, report)
        print(
            f"PASS m_max=300 rss={peak_rss:.2f}GiB "
            f"max_sample={report['max_sample_seconds']:.1f}s gate={gate_path}"
        )
    else:
        gate_path.unlink(missing_ok=True)
        print("FAIL " + "; ".join(reasons))
        raise SystemExit(1)


if __name__ == "__main__":
    main()
