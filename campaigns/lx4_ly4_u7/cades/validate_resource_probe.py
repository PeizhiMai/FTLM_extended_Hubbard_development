#!/usr/bin/env python3
"""Apply the 4x4 memory/runtime gate to twist 005's largest block."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path


DEFAULT_ROOT = Path("/lustre/or-scratch24/scratch/9pm/ftlm_codex/lx4_ly4_u7_campaign")
EXPECTED_DIMENSION = 10_353_252
EXPECTED_SAMPLES = set(range(16))
MAX_SAMPLE_SECONDS = 180 * 60


def atomic_json(path: Path, value):
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--campaign-root", type=Path, default=DEFAULT_ROOT)
    args = parser.parse_args()
    root = args.campaign_root.resolve()
    progress = root / "runs/twist_005/progress.jsonl"
    if not progress.exists():
        raise SystemExit(f"missing progress log: {progress}")

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
        peak_rss = max(peak_rss, float(event.get("max_rss_gb", 0.0)))
        block = event.get("current_block")
        if block and (
            block.get("n_up"), block.get("n_down"), block.get("mx"), block.get("my")
        ) == (8, 8, 0, 0):
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
        reasons.append(f"missing durable timing for sample IDs {missing_samples}")
    slow = {sample: seconds for sample, seconds in sample_durations.items() if seconds >= MAX_SAMPLE_SECONDS}
    if slow:
        reasons.append(f"samples exceeded 180 minutes: {slow}")
    if peak_rss > 210.0:
        reasons.append(f"peak RSS {peak_rss:.2f} GiB exceeds 210 GiB stop threshold")

    production_threads = 8 if 180.0 <= peak_rss <= 210.0 else 16
    report = {
        "passed": not reasons,
        "expected_dimension": EXPECTED_DIMENSION,
        "observed_dimensions": sorted(dimensions),
        "peak_rss_gb": peak_rss,
        "sample_durations_seconds": sample_durations,
        "max_sample_seconds": max(sample_durations.values(), default=None),
        "missing_sample_ids": missing_samples,
        "production_threads": production_threads,
        "memory_mode": "batch8" if production_threads == 8 else "samples16",
        "malformed_progress_lines": malformed,
        "reasons": reasons,
    }
    report_path = root / "resource_probe_report.json"
    atomic_json(report_path, report)
    gate_path = root / "RESOURCE_GATE_PASSED.json"
    if report["passed"]:
        atomic_json(gate_path, report)
        print(
            f"PASS rss={peak_rss:.2f}GiB max_sample={report['max_sample_seconds']:.1f}s "
            f"production_threads={production_threads} gate={gate_path}"
        )
    else:
        gate_path.unlink(missing_ok=True)
        print("FAIL " + "; ".join(reasons))
        raise SystemExit(1)


if __name__ == "__main__":
    main()
