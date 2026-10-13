#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$project_root/scripts/r_environment.sh"

if [[ "${1:-}" != "Rscript" ]]; then
  echo "Usage: scripts/run_r.sh Rscript <script> [args...]" >&2
  exit 2
fi
shift

if resolve_r_environment; then
  :
else
  status=$?
  if [[ "$status" == 3 ]]; then
    echo "Missing '$environment_name' R environment. Run: make setup-r" >&2
  fi
  exit 1
fi

rscript="$environment_prefix/bin/Rscript"
if [[ ! -x "$rscript" ]]; then
  echo "Missing '$environment_name' R environment. Run: make setup-r" >&2
  exit 1
fi

cd "$project_root"
isolate_r_environment
exec "$rscript" "$@"
