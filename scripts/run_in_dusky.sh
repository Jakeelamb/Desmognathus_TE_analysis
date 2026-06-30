#!/usr/bin/env bash
# Run a command inside the repo's Dusky conda environment.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
ENV_NAME="${DUSKY_CONDA_ENV:-Dusky}"

configure_env() {
    export PATH="$CONDA_PREFIX/bin:$PATH"
    export DUSKY_ENV_WRAPPED=1

    # Keep R scripts reproducible: user startup files can inject host-level
    # libraries built for a different R version than the conda environment.
    export R_ENVIRON_USER=/dev/null
    export R_PROFILE_USER=/dev/null
    export R_LIBS_USER="$CONDA_PREFIX/lib/R/library"
    export R_LIBS_SITE="$CONDA_PREFIX/lib/R/library"
}

if [ "$#" -eq 0 ]; then
    cat >&2 <<EOF
Usage: scripts/run_in_dusky.sh <command> [args...]

Examples:
  scripts/run_in_dusky.sh python verify_setup.py --skip-data
  scripts/run_in_dusky.sh python scripts/processing/dnaPipe.py
  scripts/run_in_dusky.sh Rscript path_analysis/scripts/path_model_scaffold.R --summary-only
EOF
    exit 2
fi

if [ -n "${CONDA_PREFIX:-}" ] && [ "$(basename "$CONDA_PREFIX")" = "$ENV_NAME" ]; then
    configure_env
    cd "$PROJECT_ROOT"
    exec "$@"
fi

if ! command -v conda >/dev/null 2>&1; then
    echo "Error: conda is required to run the $ENV_NAME environment." >&2
    echo "Install Miniconda or source your conda initialization script first." >&2
    exit 1
fi

exec conda run -n "$ENV_NAME" bash -lc \
    'export PATH="$CONDA_PREFIX/bin:$PATH";
     export DUSKY_ENV_WRAPPED=1;
     export R_ENVIRON_USER=/dev/null;
     export R_PROFILE_USER=/dev/null;
     export R_LIBS_USER="$CONDA_PREFIX/lib/R/library";
     export R_LIBS_SITE="$CONDA_PREFIX/lib/R/library";
     cd "$1";
     shift;
     exec "$@"' \
    bash "$PROJECT_ROOT" "$@"
