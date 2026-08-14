#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT=$(cd "$(dirname "$0")/../../../.." && pwd)
REMOTE=${CADES_HOST:-cades}
CAMPAIGN_ROOT=${CAMPAIGN_ROOT:-/lustre/or-scratch24/scratch/9pm/ftlm_codex/lx4_ly4_u7_campaign}
TOOLCHAIN=${1:-intel}

cd "$REPO_ROOT"
if [[ -n $(git status --porcelain) ]]; then
  echo "Refusing to deploy a dirty tree; commit the campaign implementation first." >&2
  exit 2
fi
COMMIT=$(git rev-parse HEAD)
REMOTE_SOURCE="$CAMPAIGN_ROOT/sources/$COMMIT"
REMOTE_BIN="$CAMPAIGN_ROOT/bin-$COMMIT-$TOOLCHAIN"

ssh "$REMOTE" "mkdir -p '$REMOTE_SOURCE' '$REMOTE_BIN' '$CAMPAIGN_ROOT/runs' '$CAMPAIGN_ROOT/logs'"
git archive --format=tar HEAD | ssh "$REMOTE" "tar -xf - -C '$REMOTE_SOURCE'"
printf '%s\n' "$COMMIT" | ssh "$REMOTE" "cat > '$REMOTE_SOURCE/GIT_COMMIT'"
ssh "$REMOTE" \
  "bash '$REMOTE_SOURCE/campaigns/lx4_ly4_u7/cades/build_cades.sh' '$REMOTE_SOURCE' '$REMOTE_BIN' '$TOOLCHAIN'"
ssh "$REMOTE" "ln -sfn '$REMOTE_BIN' '$CAMPAIGN_ROOT/bin'; ln -sfn '$REMOTE_SOURCE' '$CAMPAIGN_ROOT/source'"
echo "deployed_commit=$COMMIT"
echo "campaign_root=$CAMPAIGN_ROOT"
echo "bin=$REMOTE_BIN"
