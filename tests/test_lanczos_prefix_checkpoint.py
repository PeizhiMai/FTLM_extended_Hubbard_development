#!/usr/bin/env python3
"""Regression tests for one-recurrence, multi-prefix V3 checkpoints."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import struct
import subprocess
import tempfile
from pathlib import Path


V3_MAGIC = b"FTLMCP3\n"
RECORD_MAGIC = b"F3REC1\n\0"
HEADER = struct.Struct("=IQQ")


def run(command: list[str], expected: int = 0) -> subprocess.CompletedProcess:
    result = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != expected:
        raise AssertionError(
            f"return code {result.returncode}, expected {expected}: {' '.join(command)}\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    return result


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def records(path: Path):
    raw = path.read_bytes()
    assert raw.startswith(V3_MAGIC)
    cursor = len(V3_MAGIC)
    found = []
    while cursor < len(raw):
        assert raw[cursor : cursor + 8] == RECORD_MAGIC
        kind, size, checksum = HEADER.unpack_from(raw, cursor + 8)
        begin = cursor + 8 + HEADER.size
        end = begin + size
        assert end <= len(raw)
        found.append((kind, raw[begin:end], checksum))
        cursor = end
    return found


def sample_depths(path: Path):
    found: dict[tuple[int, int, int, int, int], set[int]] = {}
    for kind, payload, _ in records(path):
        if kind != 1:
            continue
        n_up, n_down, mx, my = struct.unpack_from("=iiii", payload, 0)
        sample_id = struct.unpack_from("=i", payload, 28)[0]
        lanczos_steps = struct.unpack_from("=i", payload, 32)[0]
        found.setdefault((n_up, n_down, mx, my, sample_id), set()).add(lanczos_steps)
    return found


def command(
    ftlm: Path,
    checkpoint: Path,
    output: Path,
    samples: int,
    max_steps: int,
    save_steps: str | None,
):
    result = [
        str(ftlm),
        "--lx", "2", "--ly", "2",
        "--u", "7", "--tp", "0",
        "--phix", "0.25", "--phiy", "0.25",
        "--beta-list", "2.857142857142857,12.5",
        "--mu-min", "-3", "--mu-max", "4", "--mu-count", "21",
        "--lanczos-steps", str(max_steps),
        "--exact-block-threshold", "0",
        "--threads", "2", "--seed", "12345",
        "--checkpoint", str(checkpoint),
        "--samples", str(samples),
        "--output", str(output),
    ]
    if save_steps is not None:
        result.extend(["--lanczos-save-steps", save_steps])
    return result


def reduce(reducer: Path, checkpoint: Path, output: Path, samples: int, steps: int):
    run([
        str(reducer),
        "--checkpoint", str(checkpoint),
        "--samples", str(samples),
        "--lanczos-steps", str(steps),
        "--beta-list", "2.857142857142857,12.5",
        "--mu-min", "-3", "--mu-max", "4", "--mu-count", "21",
        "--output", str(output),
    ])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ftlm", required=True, type=Path)
    parser.add_argument("--reducer", required=True, type=Path)
    parser.add_argument("--keep", action="store_true")
    args = parser.parse_args()
    root = Path(tempfile.mkdtemp(prefix="ftlm-prefix-test-"))
    try:
        checkpoint = root / "prefix.ftlmcp"
        direct12 = root / "direct-m12-r2.csv"
        first = run(command(args.ftlm, checkpoint, direct12, 2, 12, "4,8,12"))
        assert "new_samples=168" in first.stdout
        depths = sample_depths(checkpoint)
        assert len(depths) == 168
        assert all(value == {4, 8, 12} for value in depths.values())

        # Every prefix must equal a standalone run at that depth: the same
        # random vector and recurrence are used, with no independent sampling.
        prefix_hashes = {}
        for steps in (4, 8, 12):
            reduced = root / f"prefix-m{steps}-r2.csv"
            reduce(args.reducer, checkpoint, reduced, 2, steps)
            prefix_hashes[steps] = sha256(reduced)
            fresh_cp = root / f"fresh-m{steps}.ftlmcp"
            fresh = root / f"fresh-m{steps}.csv"
            run(command(args.ftlm, fresh_cp, fresh, 2, steps, None))
            assert sha256(reduced) == sha256(fresh)
        assert prefix_hashes[12] == sha256(direct12)

        # Batch mode reads the checkpoint once and emits every requested
        # (m,R) reduction, which is the production job-wrapper path.
        run([
            str(args.reducer), "--checkpoint", str(checkpoint),
            "--samples-list", "1,2", "--lanczos-steps-list", "4,8,12",
            "--beta-list", "2.857142857142857,12.5",
            "--mu-min", "-3", "--mu-max", "4", "--mu-count", "21",
            "--output-template", str(root / "batch-m{m}-R{R}.csv"),
        ])
        for steps, expected_hash in prefix_hashes.items():
            assert sha256(root / f"batch-m{steps:03d}-R002.csv") == expected_hash
        assert all((root / f"batch-m{steps:03d}-R001.csv").exists() for steps in (4, 8, 12))

        # Increasing R appends only the permanent sample ID 2 at all prefixes.
        extended = root / "direct-m12-r3.csv"
        extension = run(command(args.ftlm, checkpoint, extended, 3, 12, "4,8,12"))
        assert "reused_samples=168" in extension.stdout
        assert "new_samples=84" in extension.stdout
        depths = sample_depths(checkpoint)
        assert len(depths) == 252
        assert all(value == {4, 8, 12} for value in depths.values())
        for steps, old_hash in prefix_hashes.items():
            after = root / f"prefix-m{steps}-r2-after-r3.csv"
            reduce(args.reducer, checkpoint, after, 2, steps)
            assert sha256(after) == old_hash

        # V3 deliberately leaves m_max out of immutable metadata.  A later
        # deeper run replays each recurrence but appends only the missing m=16
        # record; all earlier reductions remain bit-for-bit stable.
        direct16 = root / "direct-m16-r3.csv"
        run(command(args.ftlm, checkpoint, direct16, 3, 16, "4,8,12,16"))
        assert all(value == {4, 8, 12, 16} for value in sample_depths(checkpoint).values())
        for steps, old_hash in prefix_hashes.items():
            after = root / f"prefix-m{steps}-r2-after-m16.csv"
            reduce(args.reducer, checkpoint, after, 2, steps)
            assert sha256(after) == old_hash

        fresh16_cp = root / "fresh-m16.ftlmcp"
        fresh16 = root / "fresh-m16.csv"
        run(command(args.ftlm, fresh16_cp, fresh16, 3, 16, None))
        assert sha256(fresh16) == sha256(direct16)

        missing_depth = run([
            str(args.reducer), "--checkpoint", str(checkpoint), "--samples", "3", "--status"
        ], expected=1)
        assert "--lanczos-steps M is required" in missing_depth.stderr
        status_path = root / "status.json"
        run([
            str(args.reducer), "--checkpoint", str(checkpoint), "--samples", "3",
            "--lanczos-steps", "16", "--status-json", str(status_path),
        ])
        status = json.loads(status_path.read_text())
        assert status["complete"] and status["target_m"] == 16
        print(f"multi-prefix checkpoint tests passed: {root}")
    finally:
        if args.keep:
            print(f"kept test directory: {root}")
        else:
            shutil.rmtree(root)


if __name__ == "__main__":
    main()
