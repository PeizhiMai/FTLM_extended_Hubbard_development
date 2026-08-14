#!/usr/bin/env bash
#SBATCH --account=cnms
#SBATCH --partition=high_mem
#SBATCH --qos=cnms-high_mem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=16G
#SBATCH --time=00:30:00

set -eo pipefail
source /etc/profile
set -u
module purge
module load intel/2024.1.0 mkl/2024.1.0

: "${CAMPAIGN_ROOT:?}"
BIN="$CAMPAIGN_ROOT/bin"
OUT="$CAMPAIGN_ROOT/validation"
mkdir -p "$OUT"
export OMP_PROC_BIND=close OMP_PLACES=cores OMP_DYNAMIC=false
export MKL_NUM_THREADS=1 MKL_DYNAMIC=false

COMMON_3X3=(
  --lx 3 --ly 3 --tx 1 --ty 1 --tp 0 --u 7 --v 0
  --phix .25 --phiy .25
  --beta-list 2.857142857142857,12.5
  --mu-min -3 --mu-max 4 --mu-count 281
  --lanczos-steps 80 --exact-block-threshold 256
  --threads 32 --seed 12345 --max-sector-dim 1000000
  --checkpoint "$OUT/lx3_ly3_v2.ftlmcp"
)
"$BIN/ftlm_n_vs_mu" "${COMMON_3X3[@]}" --samples 16 --output "$OUT/lx3_ly3_R016_before.csv"
sha256sum "$OUT/lx3_ly3_R016_before.csv" > "$OUT/lx3_ly3_R016_before.sha256"
"$BIN/ftlm_n_vs_mu" "${COMMON_3X3[@]}" --samples 128 --output "$OUT/lx3_ly3_R128.csv"
"$BIN/ftlm_reduce_checkpoint" \
  --checkpoint "$OUT/lx3_ly3_v2.ftlmcp" --samples 16 \
  --beta-list 2.857142857142857,12.5 \
  --mu-min -3 --mu-max 4 --mu-count 281 \
  --output "$OUT/lx3_ly3_R016_after.csv"
cmp "$OUT/lx3_ly3_R016_before.csv" "$OUT/lx3_ly3_R016_after.csv"

"$BIN/ftlm_n_vs_mu" \
  --lx 3 --ly 4 --tx 1 --ty 1 --tp 0 --u 7 --v 0 \
  --phix .25 --phiy .25 \
  --beta-list 2.857142857142857,12.5 \
  --mu-min -3 --mu-max 4 --mu-count 281 \
  --samples 4 --lanczos-steps 80 --exact-block-threshold 256 \
  --threads 4 --seed 12345 --max-sector-dim 10000000 \
  --checkpoint "$OUT/lx3_ly4_v2.ftlmcp" \
  --output "$OUT/lx3_ly4_R004.csv"

python3 - "$OUT/lx3_ly4_R004.csv" <<'PY'
import csv, math, sys
from collections import defaultdict
groups=defaultdict(list)
with open(sys.argv[1],newline='') as f:
    for row in csv.DictReader(f): groups[float(row['beta'])].append(row)
assert len(groups)==2 and all(len(rows)==281 for rows in groups.values())
for beta,rows in groups.items():
    rows.sort(key=lambda r:float(r['mu']))
    n=[float(r['n']) for r in rows]; k=[float(r['compressibility']) for r in rows]
    assert all(math.isfinite(x) for x in n+k)
    assert all(b+1e-9>=a for a,b in zip(n,n[1:]))
    assert min(k)>=-1e-10
    x=[1-v for v in n]
    assert min(x)<=0.002 and max(x)>=0.35
print('3x4 v2 validation passed')
PY

sha256sum "$OUT"/*.ftlmcp "$OUT"/*.meta "$OUT"/*.csv > "$OUT/SHA256SUMS.txt"
echo "FTLM_VALIDATION_COMPLETE utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
