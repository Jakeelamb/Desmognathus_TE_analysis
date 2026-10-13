#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/activate_pixi.sh"

if [[ "${1:-}" != Rscript ]]; then
  echo "Usage: scripts/run_r.sh Rscript <script> [args...]" >&2
  exit 2
fi
shift
resolve_pixi_environment || exit 1
cd "$project_root"
if [[ "${PIXI_PROJECT_ROOT:-}" != "$project_root" || "${CONDA_PREFIX:-}" != "$environment_prefix" ]]; then
  command -v pixi >/dev/null || { echo 'Pixi is required; install Pixi and run make setup.' >&2; exit 1; }
  isolate_pixi_environment
  exec pixi run --manifest-path "$project_root/pixi.toml" --environment "$environment_name" \
    --locked --no-install --executable bash "$project_root/scripts/run_r.sh" Rscript "$@"
fi
isolate_pixi_environment
exec "$environment_prefix/bin/Rscript" "$@"
