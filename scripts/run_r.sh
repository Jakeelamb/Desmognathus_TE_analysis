#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
environment_name="${DESMOGNATHUS_R_ENV:-desmognathus}"

if [[ "${1:-}" != "Rscript" ]]; then
  echo "Usage: scripts/run_r.sh Rscript <script> [args...]" >&2
  exit 2
fi
shift

if ! command -v conda >/dev/null 2>&1; then
  echo "Conda is required for the single '$environment_name' R environment." >&2
  exit 1
fi

if [[ -n "${CONDA_PREFIX:-}" && "$(basename "$CONDA_PREFIX")" == "$environment_name" ]]; then
  environment_prefix="$CONDA_PREFIX"
else
  conda_base="$(conda info --base)"
  environment_prefix="$conda_base/envs/$environment_name"
fi

rscript="$environment_prefix/bin/Rscript"
if [[ ! -x "$rscript" ]]; then
  echo "Missing '$environment_name' R environment. Run: make setup-r" >&2
  exit 1
fi

cd "$project_root"
export CONDA_PREFIX="$environment_prefix"
export PATH="$environment_prefix/bin:$PATH"
export R_ENVIRON_USER=/dev/null
export R_PROFILE_USER=/dev/null
exec "$rscript" "$@"
