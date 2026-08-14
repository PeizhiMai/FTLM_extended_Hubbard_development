#!/usr/bin/env python3
"""Generate the fixed 4x4 twist grid and permanent per-twist seeds."""

import csv
from pathlib import Path


def main():
    output = Path(__file__).with_name("twists.csv")
    values = (0.0, 0.25, 0.5, 0.75)
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=["twist_id", "phix", "phiy", "seed"],
            lineterminator="\n",
        )
        writer.writeheader()
        for iy, phiy in enumerate(values):
            for ix, phix in enumerate(values):
                twist_id = iy * len(values) + ix
                writer.writerow(
                    {
                        "twist_id": f"{twist_id:03d}",
                        "phix": f"{phix:.2f}",
                        "phiy": f"{phiy:.2f}",
                        "seed": 12345 + 1000003 * twist_id,
                    }
                )
    print(output.resolve())


if __name__ == "__main__":
    main()
