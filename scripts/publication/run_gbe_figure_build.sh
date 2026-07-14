#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
CONDA_BASE="$(conda info --base)"
RSCRIPT="${CONDA_BASE}/envs/desmo-gbe-figures/bin/Rscript"

if [[ ! -x "${RSCRIPT}" ]]; then
  echo "Missing desmo-gbe-figures environment. Run: conda env create -f Publication/figure_environment.yml" >&2
  exit 1
fi

cd "${ROOT}"
exec "${RSCRIPT}" --no-environ scripts/publication/build_gbe_figures.R
