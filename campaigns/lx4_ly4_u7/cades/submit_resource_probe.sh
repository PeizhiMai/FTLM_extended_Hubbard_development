#!/usr/bin/env bash
set -euo pipefail

CAMPAIGN_ROOT=${CAMPAIGN_ROOT:-/lustre/or-scratch24/scratch/9pm/ftlm_codex/lx4_ly4_u7_campaign}
ACCOUNT=${ACCOUNT:-cnms}
JOB_SCRIPT="$CAMPAIGN_ROOT/source/campaigns/lx4_ly4_u7/cades/job_one_twist.sbatch"
mkdir -p "$CAMPAIGN_ROOT/logs" "$CAMPAIGN_ROOT/runs/twist_005"

active=$(squeue -h -n f4x4_probe005 -o '%i|%T|%R' | head -1 || true)
if [[ -n $active ]]; then
  echo "probe already active: $active"
  exit 0
fi

exports="ALL,CAMPAIGN_ROOT=$CAMPAIGN_ROOT,TWIST_ID=005,PHIX=0.25,PHIY=0.25,SEED=5012360,TARGET_R=16,FTLM_THREADS=16,ONLY_BLOCK_NUP=8,ONLY_BLOCK_NDOWN=8,ONLY_BLOCK_MX=0,ONLY_BLOCK_MY=0"
job_id=$(sbatch --parsable \
  --account "$ACCOUNT" \
  --job-name f4x4_probe005 \
  --output "$CAMPAIGN_ROOT/logs/%x-%j.out" \
  --error "$CAMPAIGN_ROOT/logs/%x-%j.err" \
  --export "$exports" \
  "$JOB_SCRIPT")
echo "submitted resource probe job=${job_id%%;*} checkpoint=$CAMPAIGN_ROOT/runs/twist_005/twist_005.ftlmcp"
