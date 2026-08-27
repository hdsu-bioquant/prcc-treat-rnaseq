#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat >&2 <<'USAGE'
Usage:
  tests/maintainers/update_validation_baseline.sh synthetic [--apply]
  tests/maintainers/update_validation_baseline.sh realistic RUN_DIR [--apply]

Without --apply, show the difference between the generated validation manifest and
its maintained baseline and make no repository changes.

With --apply, replace only the corresponding maintained expected validation
manifest after showing the difference. This script never updates synthetic input
fixture checksums or the canonical GDC installed-reference manifest.
USAGE
    exit 2
}

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
MODE="${1:-}"
[[ -n "$MODE" ]] || usage
shift

APPLY=false
RUN_DIR=""

case "$MODE" in
    synthetic)
        if [[ "${1:-}" == "--apply" ]]; then
            APPLY=true
            shift
        fi
        [[ $# -eq 0 ]] || usage
        CANDIDATE="$REPO_ROOT/tests/synthetic/output/results/run/validation_checksums.sha256"
        EXPECTED="$REPO_ROOT/tests/synthetic/expected/validation_checksums.sha256"
        NEXT_COMMAND="bash tests/synthetic/run_test.sh"
        ;;

    realistic)
        [[ $# -ge 1 ]] || usage
        RUN_DIR="$1"
        shift
        if [[ "${1:-}" == "--apply" ]]; then
            APPLY=true
            shift
        fi
        [[ $# -eq 0 ]] || usage
        CANDIDATE="$RUN_DIR/output/results/run/validation_checksums.sha256"
        EXPECTED="$REPO_ROOT/tests/real/expected/validation_checksums.sha256"
        NEXT_COMMAND="python tests/real/validate_results.py --run-dir \"$RUN_DIR\""
        ;;

    *)
        usage
        ;;
esac

if [[ ! -s "$CANDIDATE" ]]; then
    echo "ERROR: generated validation manifest is missing or empty:" >&2
    echo "  $CANDIDATE" >&2
    echo "Run the corresponding qualification with --skip-frozen-baseline first." >&2
    exit 1
fi

if [[ ! -s "$EXPECTED" ]]; then
    echo "ERROR: maintained validation baseline is missing or empty:" >&2
    echo "  $EXPECTED" >&2
    exit 1
fi

# Reject malformed checksum lines before presenting or copying a candidate.
validate_manifest_format() {
    local path="$1"
    awk '
        NF < 2 || $1 !~ /^[0-9a-fA-F]{64}$/ { bad=1 }
        END { exit bad ? 1 : 0 }
    ' "$path"
}

if ! validate_manifest_format "$CANDIDATE"; then
    echo "ERROR: generated manifest has malformed SHA256 entries: $CANDIDATE" >&2
    exit 1
fi
if ! validate_manifest_format "$EXPECTED"; then
    echo "ERROR: maintained baseline has malformed SHA256 entries: $EXPECTED" >&2
    exit 1
fi

echo "Generated manifest:  $CANDIDATE"
echo "Maintained baseline:  $EXPECTED"
echo

if cmp -s "$EXPECTED" "$CANDIDATE"; then
    echo "PASS: generated validation manifest already matches the maintained baseline."
    echo "No baseline update is required."
    exit 0
fi

echo "Candidate baseline differences:"
echo "------------------------------------------------------------"
diff -u "$EXPECTED" "$CANDIDATE" || true
echo "------------------------------------------------------------"
echo

if [[ "$APPLY" != true ]]; then
    echo "Preview only: no files changed."
    echo "After the differences are understood and the required reproducibility checks are complete,"
    if [[ "$MODE" == "synthetic" ]]; then
        echo "apply with: bash tests/maintainers/update_validation_baseline.sh synthetic --apply"
    else
        printf 'apply with: bash tests/maintainers/update_validation_baseline.sh realistic %q --apply\n' "$RUN_DIR"
    fi
    exit 0
fi

TMP="${EXPECTED}.tmp.$$"
trap 'rm -f "$TMP"' EXIT
cp "$CANDIDATE" "$TMP"
mv "$TMP" "$EXPECTED"
trap - EXIT

echo "UPDATED: ${EXPECTED#$REPO_ROOT/}"
echo

echo "Re-run ordinary validation against the maintained baseline:"
echo "  cd \"$REPO_ROOT\""
echo "  $NEXT_COMMAND"
