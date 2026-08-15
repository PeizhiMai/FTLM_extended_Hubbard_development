#!/usr/bin/env python3
"""Generate the D4-irreducible 8x8 midpoint grid in the 4x4 rBZ."""

from __future__ import annotations

import csv
import json
from fractions import Fraction
from pathlib import Path


LATTICE_LENGTH = 4
GRID_SIZE = 8
FIRST_TWIST_ID = 800
BASE_SEED = 12345
SEED_STRIDE = 1000003


def positive_midpoints() -> list[Fraction]:
    """Positive half of the centered midpoint grid in boundary-flux units."""

    return [Fraction(2 * index + 1, 2 * GRID_SIZE) for index in range(GRID_SIZE // 2)]


def main():
    directory = Path(__file__).resolve().parent
    output = directory / "twists.csv"
    metadata_output = directory / "twist_quadrature.json"
    levels = positive_midpoints()
    rows = []
    twist_index = 0
    for iy, phiy in enumerate(levels):
        for ix, phix in enumerate(levels[: iy + 1]):
            twist_id = FIRST_TWIST_ID + twist_index
            # The executable uses k_tilde = 2*pi*phi/L. D4 produces four
            # sign-related diagonal points or eight sign/exchange-related
            # off-diagonal points from each representative in this open wedge.
            multiplicity = 4 if ix == iy else 8
            rows.append(
                {
                    "twist_id": f"{twist_id:03d}",
                    "phix": f"{float(phix):.8f}",
                    "phiy": f"{float(phiy):.8f}",
                    "kx_over_pi": f"{float(2 * phix / LATTICE_LENGTH):.8f}",
                    "ky_over_pi": f"{float(2 * phiy / LATTICE_LENGTH):.8f}",
                    "seed": BASE_SEED + SEED_STRIDE * twist_id,
                    "multiplicity": multiplicity,
                }
            )
            twist_index += 1

    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    metadata = {
        "schema_version": 1,
        "tag": "rbz8_midpoint",
        "description": "8x8 midpoint quadrature of the 4x4-cluster reduced Brillouin zone",
        "cluster_lx": LATTICE_LENGTH,
        "cluster_ly": LATTICE_LENGTH,
        "grid_nx": GRID_SIZE,
        "grid_ny": GRID_SIZE,
        "full_domain": "-pi/4 <= kx,ky < pi/4",
        "representative_domain": "0 < kx <= ky < pi/4",
        "symmetry": "D4: independent sign reflections and x-y exchange",
        "input_convention": "kx=2*pi*phix/lx; ky=2*pi*phiy/ly",
        "representatives": len(rows),
        "effective_twists": sum(int(row["multiplicity"]) for row in rows),
        "midpoint_grid": True,
    }
    metadata_output.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    print(output)
    print(metadata_output)


if __name__ == "__main__":
    main()
