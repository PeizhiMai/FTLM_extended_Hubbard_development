#!/usr/bin/env bash
set -euo pipefail

CAMPAIGN_ROOT=${CAMPAIGN_ROOT:-/lustre/or-scratch24/scratch/9pm/ftlm_codex/lx4_ly4_u7_campaign}
ACCOUNT=${ACCOUNT:-cnms}
SUBMIT="$CAMPAIGN_ROOT/source/campaigns/lx4_ly4_u7/cades/submit_wave.py"

# This is deliberately a separate Slurm job and a separate V3 checkpoint from
# the active legacy m=80 pilot.  The 16 largest-block samples become reusable
# records when twist 005 is subsequently continued through all blocks.
/usr/bin/python3.11 "$SUBMIT" \
  --campaign-root "$CAMPAIGN_ROOT" \
  --account "$ACCOUNT" \
  --samples 16 \
  --twist-id 005 \
  --checkpoint-series mext \
  --lanczos-max-steps 300 \
  --lanczos-save-steps 80,120,160,200,250,300 \
  --only-block 8,8,0,0 \
  --allow-ungated
