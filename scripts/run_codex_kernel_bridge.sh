#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

source /software/miniconda3/bin/activate base
conda activate pyscenic

exec python -u "${project_root}/scripts/codex_kernel_bridge.py" \
  --project-root "${project_root}" \
  "$@"
