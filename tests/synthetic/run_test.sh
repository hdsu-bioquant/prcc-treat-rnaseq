#!/usr/bin/env bash

set -euo pipefail

# Always operate from repository root
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

echo "Cleaning previous synthetic test outputs..."
rm -rf tests/synthetic/results
rm -rf tests/synthetic/reference/star_index

echo "Running synthetic Snakemake smoke test..."
snakemake \
    --snakefile workflow/Snakefile \
    --configfile tests/synthetic/config.yaml \
    --cores 2 \
    --software-deployment-method apptainer \
    --printshellcmds

echo "Validating synthetic results..."
python tests/synthetic/validate_results.py