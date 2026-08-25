#!/usr/bin/env bash
# Download/verify the exact compressed ENA FASTQs used for realistic qualification.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TEST_DIR="$REPO_ROOT/tests/real"
MANIFEST="$TEST_DIR/data_manifest.tsv"
DATA_DIR="$TEST_DIR/data"
VERIFY_ONLY=false

usage() {
    cat <<EOF_USAGE
Usage: bash tests/real/get_test_data.sh [--verify-only]

Without options, download missing/invalid files using the pinned URLs in
$MANIFEST, then verify exact byte size, published MD5 and gzip integrity.

--verify-only  Never access the network; require all three files to be present
               and valid.
EOF_USAGE
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --verify-only) VERIFY_ONLY=true ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
    shift
done

for cmd in md5sum gzip stat awk basename; do
    command -v "$cmd" >/dev/null 2>&1 || {
        echo "ERROR: required command not found: $cmd" >&2
        exit 1
    }
done
if [[ "$VERIFY_ONLY" == false ]]; then
    command -v wget >/dev/null 2>&1 || {
        echo "ERROR: wget is required for downloading qualification FASTQs" >&2
        exit 1
    }
fi

[[ -r "$MANIFEST" ]] || { echo "ERROR: manifest not readable: $MANIFEST" >&2; exit 1; }
mkdir -p "$DATA_DIR"

expected_header=$'run_accession\tlibrary_id\tread_role\tfilename\turl\tbytes\tmd5'
observed_header="$(head -n 1 "$MANIFEST")"
if [[ "$observed_header" != "$expected_header" ]]; then
    echo "ERROR: unexpected header in $MANIFEST" >&2
    echo "Observed: $observed_header" >&2
    exit 1
fi

verify_one() {
    local path="$1" expected_bytes="$2" expected_md5="$3"
    [[ -s "$path" ]] || return 1

    local observed_bytes
    observed_bytes="$(stat -c '%s' "$path")"
    [[ "$observed_bytes" == "$expected_bytes" ]] || return 1

    local observed_md5
    observed_md5="$(md5sum "$path" | awk '{print $1}')"
    [[ "$observed_md5" == "$expected_md5" ]] || return 1

    gzip -t "$path" >/dev/null 2>&1 || return 1
    return 0
}

n=0
while IFS=$'\t' read -r run_accession library_id read_role filename url bytes md5; do
    [[ "$run_accession" == "run_accession" ]] && continue
    [[ -z "$run_accession" ]] && continue
    n=$((n + 1))

    if [[ "$(basename "$url")" != "$filename" ]]; then
        echo "ERROR: manifest filename does not match pinned URL basename: $filename" >&2
        exit 1
    fi
    if [[ ! "$bytes" =~ ^[0-9]+$ || ! "$md5" =~ ^[0-9a-f]{32}$ ]]; then
        echo "ERROR: malformed byte/MD5 metadata for $filename" >&2
        exit 1
    fi

    target="$DATA_DIR/$filename"
    if verify_one "$target" "$bytes" "$md5"; then
        echo "PASS: already valid: $filename"
        continue
    fi

    if [[ "$VERIFY_ONLY" == true ]]; then
        echo "ERROR: missing or invalid qualification FASTQ: $target" >&2
        echo "Run first: bash tests/real/get_test_data.sh" >&2
        exit 1
    fi

    if [[ -e "$target" ]]; then
        observed_bytes="$(stat -c '%s' "$target" 2>/dev/null || echo 0)"
        if [[ "$observed_bytes" =~ ^[0-9]+$ ]] && (( observed_bytes > bytes )); then
            echo ">> Removing oversized/invalid local file before retry: $filename"
            rm -f "$target"
        elif [[ "$observed_bytes" == "$bytes" ]]; then
            echo ">> Removing complete-size file with failed integrity check: $filename"
            rm -f "$target"
        else
            echo ">> Resuming partial download: $filename"
        fi
    else
        echo ">> Downloading: $filename"
    fi

    download_one() {
        (
            cd "$DATA_DIR"
            wget \
                --continue \
                --tries=10 \
                --timeout=30 \
                --read-timeout=30 \
                --waitretry=5 \
                --retry-connrefused \
                "$url"
        )
    }

    download_one

    if ! verify_one "$target" "$bytes" "$md5"; then
        echo ">> Verification failed after resume/download; retrying once from scratch: $filename"
        rm -f "$target"
        download_one
    fi

    if ! verify_one "$target" "$bytes" "$md5"; then
        echo "ERROR: downloaded file failed size/MD5/gzip verification: $filename" >&2
        exit 1
    fi
    echo "PASS: downloaded and verified: $filename"
done < "$MANIFEST"

if [[ "$n" -ne 3 ]]; then
    echo "ERROR: expected exactly 3 qualification FASTQ entries; observed $n" >&2
    exit 1
fi

echo "PASS: realistic qualification FASTQ set is complete and exact: $DATA_DIR"
