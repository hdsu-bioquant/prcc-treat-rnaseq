#!/usr/bin/env bash

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

TEST_DIR="tests/synthetic"
LOG_DIR="$TEST_DIR/logs"
TIMESTAMP="$(date '+%Y%m%d_%H%M%S')"
LOG_FILE="$LOG_DIR/run_test_${TIMESTAMP}.log"

mkdir -p "$LOG_DIR"

# Show output in the terminal and save stdout + stderr to a local log.
exec > >(tee "$LOG_FILE") 2>&1

finish() {
    rc=$?
    trap - EXIT
    set +e

    echo
    echo "============================================================"

    if [[ $rc -eq 0 ]]; then
        echo "Synthetic smoke test PASSED."
    else
        echo "Synthetic smoke test FAILED (exit code: $rc)."
    fi

    echo "Log saved to: $LOG_FILE"
    echo "============================================================"

    exit "$rc"
}

trap finish EXIT

command_version() {
    local label="$1"
    shift

    if command -v "$1" >/dev/null 2>&1; then
        printf '%-14s %s\n' "$label:" "$("$@" 2>&1 | head -n 1)"
    else
        printf '%-14s %s\n' "$label:" "unavailable"
    fi
}

echo "============================================================"
echo "pRCC-RNA-Seq synthetic smoke test"
echo "============================================================"
echo "Date:          $(date --iso-8601=seconds 2>/dev/null || date)"
echo "Host:          $(hostname 2>/dev/null || echo unknown)"

command_version "Snakemake" snakemake --version
command_version "Python" python --version

if command -v apptainer >/dev/null 2>&1; then
    command_version "Apptainer" apptainer --version
elif command -v singularity >/dev/null 2>&1; then
    command_version "Singularity" singularity --version
else
    printf '%-14s %s\n' "Container:" "apptainer/singularity unavailable"
fi

echo "Log file:      $LOG_FILE"
echo "============================================================"
echo

echo "Verifying synthetic test fixture..."
sha256sum --check --quiet "$TEST_DIR/checksums.sha256"
echo "PASS: synthetic fixture checksums"
echo

echo "Cleaning previous synthetic test outputs..."
rm -rf "$TEST_DIR/results"
rm -rf "$TEST_DIR/reference/star_index"
echo

echo "Running synthetic Snakemake smoke test..."
snakemake \
    --snakefile workflow/Snakefile \
    --configfile "$TEST_DIR/config.yaml" \
    --cores 2 \
    --software-deployment-method apptainer \
    --printshellcmds

echo
echo "Validating synthetic results..."
python "$TEST_DIR/validate_results.py"