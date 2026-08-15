#!/usr/bin/env python3
"""Generate the D4-irreducible pi/40 Gamma-centered 4x4 rBZ grid."""

from __future__ import annotations

import csv
import json
from pathlib import Path


LATTICE_LENGTH = 4
GRID_SIZE = 20
FIRST_NEW_TWIST_ID = 100
BASE_SEED = 12345
SEED_STRIDE = 1000003
PRESERVED_TWISTS = {
    (0, 0): (0, 12345),
    (5, 5): (5, 5012360),
}


def d4_orbit(ix: int, iy: int) -> set[tuple[int, int]]:
    """D4 orbit on the periodic 20x20 grid, represented modulo 20."""

    points = set()
    for sx in (-1, 1):
        for sy in (-1, 1):
            points.add(((sx * ix) % GRID_SIZE, (sy * iy) % GRID_SIZE))
            points.add(((sx * iy) % GRID_SIZE, (sy * ix) % GRID_SIZE))
    return points


def main():
    directory = Path(__file__).resolve().parent
    output = directory / "twists.csv"
    metadata_output = directory / "twist_quadrature.json"
    rows = []
    next_twist_id = FIRST_NEW_TWIST_ID
    # In units of pi/40, the irreducible representatives are
    # 0 <= ix <= iy <= 10.  Index 10 is the periodic rBZ boundary: +pi/4
    # and -pi/4 are the same grid point modulo the cluster reciprocal vector.
    for iy in range(GRID_SIZE // 2 + 1):
        for ix in range(iy + 1):
            if (ix, iy) in PRESERVED_TWISTS:
                twist_id, seed = PRESERVED_TWISTS[(ix, iy)]
            else:
                twist_id = next_twist_id
                seed = BASE_SEED + SEED_STRIDE * twist_id
                next_twist_id += 1
            phix = ix / GRID_SIZE
            phiy = iy / GRID_SIZE
            multiplicity = len(d4_orbit(ix, iy))
            rows.append(
                {
                    "twist_id": f"{twist_id:03d}",
                    "grid_ix": ix,
                    "grid_iy": iy,
                    "phix": f"{phix:.8f}",
                    "phiy": f"{phiy:.8f}",
                    "kx_over_pi": f"{ix / 40.0:.8f}",
                    "ky_over_pi": f"{iy / 40.0:.8f}",
                    "seed": seed,
                    "multiplicity": multiplicity,
                }
            )

    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    metadata = {
        "schema_version": 1,
        "tag": "rbz20_gamma_dk_pi40",
        "description": "20x20 Gamma-centered quadrature with spacing pi/40 in the 4x4-cluster reduced Brillouin zone",
        "cluster_lx": LATTICE_LENGTH,
        "cluster_ly": LATTICE_LENGTH,
        "grid_nx": GRID_SIZE,
        "grid_ny": GRID_SIZE,
        "spacing_over_pi": 1 / 40,
        "full_domain": "-pi/4 <= kx,ky < pi/4",
        "representative_domain": "0 <= kx <= ky <= pi/4; +pi/4 represents the periodic -pi/4 boundary",
        "symmetry": "D4: independent sign reflections and x-y exchange",
        "input_convention": "kx=2*pi*phix/lx; ky=2*pi*phiy/ly",
        "representatives": len(rows),
        "effective_twists": sum(int(row["multiplicity"]) for row in rows),
        "gamma_centered": True,
        "includes_zero": True,
        "periodic_boundary_identification": "+pi/4 == -pi/4 modulo pi/2",
        "preserved_twist_ids": [
            {"twist_id": f"{twist_id:03d}", "grid_ix": pair[0], "grid_iy": pair[1]}
            for pair, (twist_id, _) in PRESERVED_TWISTS.items()
        ],
    }
    metadata_output.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    print(output)
    print(metadata_output)


if __name__ == "__main__":
    main()
