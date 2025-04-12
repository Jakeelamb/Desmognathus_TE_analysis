#!/bin/bash
# Run unit tests for the Desmognathus_TE project

set -e  # Exit on error

echo "Running Desmognathus_TE Unit Tests"
echo "=================================="

# Ensure we're in the project root
if [ ! -d "scripts" ] || [ ! -d "config" ]; then
    echo "Error: This script must be run from the project root directory."
    echo "Current directory: $(pwd)"
    exit 1
fi

# Check for Python
if ! command -v python &> /dev/null; then
    echo "Error: Python is required but not found."
    exit 1
fi

# Create test output directory if it doesn't exist
mkdir -p results/test_reports

# Run Python tests
echo "Running Python tests..."
PYTHONPATH="$PWD" python -m unittest discover -s scripts/python/tests

# Run R tests if available and if R is installed
if command -v Rscript &> /dev/null; then
    if [ -d "scripts/R/tests" ]; then
        echo "Running R tests..."
        Rscript -e "testthat::test_dir('scripts/R/tests')"
    else
        echo "No R tests found in scripts/R/tests."
    fi
else
    echo "Skipping R tests: R is not installed or not in PATH."
fi

echo "All tests completed!" 