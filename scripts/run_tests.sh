#!/bin/bash
# Run unit tests for the Desmognathus_TE project

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

cd "$PROJECT_ROOT"

if [ "${DUSKY_ENV_WRAPPED:-}" != "1" ]; then
    exec "$PROJECT_ROOT/scripts/run_in_dusky.sh" bash "$PROJECT_ROOT/scripts/run_tests.sh"
fi

echo "Running Desmognathus_TE Unit Tests" >&2
echo "==================================" >&2

# Check for Python
if ! command -v python &> /dev/null; then
    echo "Error: Python is required but not found." >&2
    exit 1
fi

# Run Python tests
echo "Running Python tests..." >&2
PYTHONPATH="$PWD" python -m unittest discover -s scripts/python/tests

# Run R tests if available and if R is installed
if command -v Rscript &> /dev/null; then
    if [ -d "scripts/R/tests" ]; then
        echo "Running R tests..." >&2
        Rscript -e "testthat::test_dir('scripts/R/tests')"
    else
        echo "No R tests found in scripts/R/tests." >&2
    fi
else
    echo "Skipping R tests: R is not installed or not in PATH." >&2
fi

echo "All tests completed!" >&2
