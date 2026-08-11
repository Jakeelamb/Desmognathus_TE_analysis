#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$project_root"
exec scripts/run_r.sh Rscript scripts/publication/build_gbe_figures.R "$@"
