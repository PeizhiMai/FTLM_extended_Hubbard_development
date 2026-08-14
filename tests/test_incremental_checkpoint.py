#!/usr/bin/env python3
"""End-to-end regression tests for extensible thermodynamic checkpoints."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import signal
import struct
import subprocess
import tempfile
import time
from pathlib import Path


V2_MAGIC = b"FTLMCP2\n"
RECORD_MAGIC = b"F2REC1\n\0"
HEADER = struct.Struct("=IQQ")


def run(command: list[str], *, expected: int = 0, stdout: Path | None = None) -> subprocess.CompletedProcess:
    result = subprocess.run(
        command,
        text=True,
        stdout=subprocess.PIPE if stdout is None else stdout.open("w"),
        stderr=subprocess.PIPE,
        check=False,
    )
    if result.returncode != expected:
        output = result.stdout if isinstance(result.stdout, str) else ""
        raise AssertionError(
            f"command returned {result.returncode}, expected {expected}: {' '.join(command)}\n"
            f"stdout:\n{output}\nstderr:\n{result.stderr}"
        )
    return result


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fnv1a64(payload: bytes) -> int:
    value = 1469598103934665603
    for byte in payload:
        value ^= byte
        value = (value * 1099511628211) & ((1 << 64) - 1)
    return value


def records(path: Path) -> list[tuple[int, bytes, int]]:
    raw = path.read_bytes()
    assert raw.startswith(V2_MAGIC)
    found: list[tuple[int, bytes, int]] = []
    cursor = len(V2_MAGIC)
    while cursor < len(raw):
        assert raw[cursor : cursor + 8] == RECORD_MAGIC
        kind, size, checksum = HEADER.unpack_from(raw, cursor + 8)
        begin = cursor + 8 + HEADER.size
        end = begin + size
        payload = raw[begin:end]
        assert len(payload) == size
        assert fnv1a64(payload) == checksum
        found.append((kind, payload, checksum))
        cursor = end
    return found


def sample_ids(path: Path) -> dict[tuple[int, int, int, int], set[int]]:
    result: dict[tuple[int, int, int, int], set[int]] = {}
    for kind, payload, _ in records(path):
        if kind != 1:
            continue
        n_up, n_down, mx, my = struct.unpack_from("=iiii", payload, 0)
        sample_id = struct.unpack_from("=i", payload, 28)[0]
        result.setdefault((n_up, n_down, mx, my), set()).add(sample_id)
    return result


def common(ftlm: Path, checkpoint: Path, output: Path, samples: int, *, threshold: int = 0) -> list[str]:
    return [
        str(ftlm),
        "--lx", "2", "--ly", "2",
        "--u", "7", "--tp", "0",
        "--phix", "0.25", "--phiy", "0.25",
        "--beta-list", "2.857142857142857,12.5",
        "--mu-min", "-3", "--mu-max", "4", "--mu-count", "21",
        "--lanczos-steps", "8",
        "--exact-block-threshold", str(threshold),
        "--threads", "2", "--seed", "12345",
        "--checkpoint", str(checkpoint),
        "--samples", str(samples),
        "--output", str(output),
    ]


def reduce_command(reducer: Path, checkpoint: Path, output: Path, samples: int, beta: str = "2.857142857142857,12.5") -> list[str]:
    return [
        str(reducer), "--checkpoint", str(checkpoint), "--samples", str(samples),
        "--beta-list", beta,
        "--mu-min", "-3", "--mu-max", "4", "--mu-count", "21",
        "--output", str(output),
    ]


def append_conflicting_duplicate(source: Path, destination: Path) -> None:
    shutil.copy2(source, destination)
    shutil.copy2(Path(str(source) + ".meta"), Path(str(destination) + ".meta"))
    first_sample = next((payload for kind, payload, _ in records(source) if kind == 1), None)
    assert first_sample is not None
    mutated = bytearray(first_sample)
    mutated[-1] ^= 1
    payload = bytes(mutated)
    with destination.open("ab") as stream:
        stream.write(RECORD_MAGIC)
        stream.write(HEADER.pack(1, len(payload), fnv1a64(payload)))
        stream.write(payload)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ftlm", required=True, type=Path)
    parser.add_argument("--reducer", required=True, type=Path)
    parser.add_argument("--keep", action="store_true")
    args = parser.parse_args()
    root = Path(tempfile.mkdtemp(prefix="ftlm-incremental-test-"))
    try:
        checkpoint = root / "incremental.ftlmcp"
        r2 = root / "r2.csv"
        run(common(args.ftlm, checkpoint, r2, 2))
        r2_hash = sha256(r2)
        ids_r2 = sample_ids(checkpoint)
        assert len(ids_r2) == 84
        assert all(ids == {0, 1} for ids in ids_r2.values())

        r4 = root / "r4-incremental.csv"
        extension = run(common(args.ftlm, checkpoint, r4, 4))
        assert "reused_samples=168" in extension.stdout
        assert "new_samples=168" in extension.stdout
        ids_r4 = sample_ids(checkpoint)
        assert all(ids_r4[key] - ids_r2[key] == {2, 3} for key in ids_r2)

        r2_after = root / "r2-after.csv"
        run(reduce_command(args.reducer, checkpoint, r2_after, 2))
        assert sha256(r2_after) == r2_hash

        fresh4_checkpoint = root / "fresh4.ftlmcp"
        fresh4 = root / "r4-fresh.csv"
        run(common(args.ftlm, fresh4_checkpoint, fresh4, 4))
        assert sha256(fresh4) == sha256(r4)

        r5 = root / "r5-incremental.csv"
        non_power = run(common(args.ftlm, checkpoint, r5, 5))
        assert "new_samples=84" in non_power.stdout
        assert all(ids == {0, 1, 2, 3, 4} for ids in sample_ids(checkpoint).values())
        fresh5_checkpoint = root / "fresh5.ftlmcp"
        fresh5 = root / "r5-fresh.csv"
        run(common(args.ftlm, fresh5_checkpoint, fresh5, 5))
        assert sha256(fresh5) == sha256(r5)

        incomplete = run(
            reduce_command(args.reducer, checkpoint, root / "must-not-exist.csv", 6),
            expected=1,
        )
        assert "incomplete for R=6" in incomplete.stderr

        status_json = root / "status.json"
        status = run([
            str(args.reducer), "--checkpoint", str(checkpoint), "--samples", "6",
            "--status", "--status-json", str(status_json),
        ])
        parsed_status = json.loads(status_json.read_text())
        assert parsed_status["next_missing"]["sample_id"] == 5
        assert len(parsed_status["missing_blocks"]) == 84
        assert "sample_ids=5" in status.stdout

        beta_direct = root / "beta-direct.csv"
        beta_reduced = root / "beta-reduced.csv"
        beta_command = common(args.ftlm, checkpoint, beta_direct, 5)
        beta_command[beta_command.index("--beta-list") + 1] = "1.5"
        run(beta_command)
        run(reduce_command(args.reducer, checkpoint, beta_reduced, 5, beta="1.5"))
        assert sha256(beta_direct) == sha256(beta_reduced)

        truncated = root / "truncated.ftlmcp"
        shutil.copy2(checkpoint, truncated)
        shutil.copy2(Path(str(checkpoint) + ".meta"), Path(str(truncated) + ".meta"))
        original_size = truncated.stat().st_size
        with truncated.open("ab") as stream:
            stream.write(b"F2R")
        repaired = run(common(args.ftlm, truncated, root / "repaired.csv", 5))
        assert "event=REPAIR_TRAILING_RECORD" in repaired.stdout
        assert truncated.stat().st_size == original_size

        conflicting = root / "conflicting.ftlmcp"
        append_conflicting_duplicate(checkpoint, conflicting)
        conflict = run([
            str(args.reducer), "--checkpoint", str(conflicting), "--samples", "5", "--status"
        ], expected=1)
        assert "Conflicting duplicate sample record" in conflict.stderr

        mismatch = common(args.ftlm, checkpoint, root / "mismatch.csv", 5)
        mismatch[mismatch.index("--phix") + 1] = "0.30"
        rejected = run(mismatch, expected=1)
        assert "metadata does not match" in rejected.stderr

        exact_checkpoint = root / "exact.ftlmcp"
        exact2 = root / "exact-r2.csv"
        exact4 = root / "exact-r4.csv"
        run(common(args.ftlm, exact_checkpoint, exact2, 2, threshold=256))
        exact_extension = run(common(args.ftlm, exact_checkpoint, exact4, 4, threshold=256))
        assert sha256(exact2) == sha256(exact4)
        assert not sample_ids(exact_checkpoint)
        assert "new_samples=0" in exact_extension.stdout

        # Interrupt an extension only after at least one durable record appears,
        # then verify restart computes only the remainder and matches a fresh run.
        interrupted_cp = root / "interrupted.ftlmcp"
        interrupted_out = root / "interrupted.csv"
        log_path = root / "interrupted.log"
        command = common(args.ftlm, interrupted_cp, interrupted_out, 16)
        with log_path.open("w") as log:
            process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, text=True)
            deadline = time.monotonic() + 10.0
            while time.monotonic() < deadline:
                if interrupted_cp.exists() and interrupted_cp.stat().st_size > 512:
                    os.kill(process.pid, signal.SIGUSR1)
                    break
                if process.poll() is not None:
                    break
                time.sleep(0.002)
            interrupted_rc = process.wait(timeout=30)
        assert interrupted_rc == 75, log_path.read_text()
        partial_status = root / "partial-status.json"
        run([
            str(args.reducer), "--checkpoint", str(interrupted_cp), "--samples", "16",
            "--status", "--status-json", str(partial_status),
        ])
        partial = json.loads(partial_status.read_text())
        assert 0 < partial["durable_sample_records"] < partial["expected_sample_records"]
        resumed = run(common(args.ftlm, interrupted_cp, interrupted_out, 16))
        assert f"reused_samples={partial['durable_sample_records']}" in resumed.stdout
        fresh16_cp = root / "fresh16.ftlmcp"
        fresh16 = root / "fresh16.csv"
        run(common(args.ftlm, fresh16_cp, fresh16, 16))
        assert sha256(fresh16) == sha256(interrupted_out)

        with r5.open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        assert len(rows) == 42
        betas = sorted({float(row["beta"]) for row in rows})
        assert len(betas) == 2 and abs(betas[0] - 2.857142857142857) < 1e-13
        assert betas[1] == 12.5
        print(f"incremental checkpoint tests passed: {root}")
    finally:
        if args.keep:
            print(f"kept test directory: {root}")
        else:
            shutil.rmtree(root)


if __name__ == "__main__":
    main()
