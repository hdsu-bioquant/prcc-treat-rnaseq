#!/usr/bin/env bash
#==============================================================================#
# Install the EXACT GDC reference resources used by pRCC-RNA-Seq:
#   - GRCh38.d1.vd1 genome FASTA
#   - GENCODE v36 annotation GTF
#   - GDC-built STAR 2.7.5c index (sjdbOverhang 100)
#
# Download URLs and official GDC MD5s live in resources/gdc_resources.tsv.
# Normal Snakemake execution never downloads these files; install them once at
# each site before production analysis.
#==============================================================================#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESOURCE_TABLE="$SCRIPT_DIR/gdc_resources.tsv"
VERIFY_SCRIPT="$SCRIPT_DIR/verify_gdc_references.sh"
DEFAULT_DEST="$SCRIPT_DIR/gdc"
CANONICAL_MANIFEST="$SCRIPT_DIR/gdc_installed_reference.sha256"
DEST=""

usage() {
    cat <<EOF_USAGE
Usage: $0 [REFERENCE_DIR]

Install and MD5-verify the exact GDC reference archives, then extract them.
Already installed targets and already valid archives are reused. If the
maintainer-frozen canonical installed-reference SHA256 manifest is available,
the extracted bundle is also qualified against it and stamped for consortium use.

Options:
  -h, --help  Show this help.

Default REFERENCE_DIR: $DEFAULT_DEST
EOF_USAGE
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            usage
            exit 0
            ;;
        --*)
            echo "ERROR: unknown option: $1" >&2
            usage >&2
            exit 2
            ;;
        *)
            if [[ -n "$DEST" ]]; then
                echo "ERROR: only one REFERENCE_DIR may be supplied" >&2
                usage >&2
                exit 2
            fi
            DEST="$1"
            ;;
    esac
    shift
done

DEST="${DEST:-$DEFAULT_DEST}"
mkdir -p "$DEST"
DEST="$(cd "$DEST" && pwd)"

if [[ ! -r "$RESOURCE_TABLE" ]]; then
    echo "ERROR: resource manifest not readable: $RESOURCE_TABLE" >&2
    exit 1
fi

archive_is_valid() {
    local archive_path="$1"
    local expected_md5="$2"
    [[ -s "$archive_path" ]] && echo "$expected_md5  $archive_path" | md5sum -c --status
}

ensure_archive() {
    local url="$1"
    local archive_path="$2"
    local expected_md5="$3"

    if archive_is_valid "$archive_path" "$expected_md5"; then
        echo ">> Reusing MD5-verified archive: $(basename "$archive_path")"
        return
    fi

    if [[ -e "$archive_path" ]]; then
        echo ">> Existing archive failed MD5 or is incomplete; replacing: $(basename "$archive_path")"
        rm -f "$archive_path"
    fi

    echo ">> Downloading: $(basename "$archive_path")"
    # Keep the direct wget transfer used successfully by the existing installer.
    wget -O "$archive_path" "$url"
    echo "$expected_md5  $archive_path" | md5sum -c -
}

install_genome() {
    local url="$1" archive="$2" md5="$3" installed="$4"
    ensure_archive "$url" "$DEST/$archive" "$md5"
    if [[ -s "$DEST/$installed" ]]; then
        echo ">> Genome FASTA already installed: $installed"
        return
    fi
    echo ">> Extracting genome FASTA"
    tar -xzf "$DEST/$archive" -C "$DEST"
    [[ -s "$DEST/$installed" ]] || {
        echo "ERROR: expected genome FASTA not found after extraction: $DEST/$installed" >&2
        exit 1
    }
}

install_gtf() {
    local url="$1" archive="$2" md5="$3" installed="$4"
    ensure_archive "$url" "$DEST/$archive" "$md5"
    if [[ -s "$DEST/$installed" ]]; then
        echo ">> Annotation GTF already installed: $installed"
        return
    fi
    echo ">> Extracting annotation GTF"
    local tmp="$DEST/.${installed}.tmp.$$"
    rm -f "$tmp"
    gzip -dc "$DEST/$archive" > "$tmp"
    mv "$tmp" "$DEST/$installed"
    [[ -s "$DEST/$installed" ]] || {
        echo "ERROR: expected annotation GTF not found after extraction: $DEST/$installed" >&2
        exit 1
    }
}

star_index_looks_complete() {
    local d="$1"
    [[ -s "$d/Genome" && -s "$d/SA" && -s "$d/SAindex" && -s "$d/genomeParameters.txt" ]]
}

install_star_index() {
    local url="$1" archive="$2" md5="$3" installed="$4"
    local target="$DEST/$installed"
    ensure_archive "$url" "$DEST/$archive" "$md5"
    if star_index_looks_complete "$target"; then
        echo ">> GDC STAR index already installed: $installed"
        return
    fi

    if [[ -e "$target" ]]; then
        echo ">> Removing incomplete STAR index before extraction: $target"
        rm -rf "$target"
    fi
    mkdir -p "$target"
    echo ">> Extracting GDC STAR index (large; this may take some time)"
    tar -xzf "$DEST/$archive" -C "$target" --strip-components=1
}

while IFS=$'\t' read -r resource_id role url archive official_md5 installed_path; do
    [[ "$resource_id" == "resource_id" ]] && continue
    [[ -z "$resource_id" ]] && continue
    case "$role" in
        genome_fasta)
            install_genome "$url" "$archive" "$official_md5" "$installed_path"
            ;;
        annotation_gtf)
            install_gtf "$url" "$archive" "$official_md5" "$installed_path"
            ;;
        star_index)
            install_star_index "$url" "$archive" "$official_md5" "$installed_path"
            ;;
        *)
            echo "ERROR: unsupported role '$role' in $RESOURCE_TABLE" >&2
            exit 1
            ;;
    esac
done < "$RESOURCE_TABLE"

echo
echo ">> Verifying installed reference structure"
bash "$VERIFY_SCRIPT" "$DEST"

if [[ -s "$CANONICAL_MANIFEST" ]]; then
    echo
    echo ">> Verifying against maintainer-frozen installed-reference SHA256 manifest"
    bash "$VERIFY_SCRIPT" --qualify "$DEST"
else
    echo
    echo ">> Canonical installed-reference SHA256 manifest is not frozen yet."
    echo ">> Skipping byte-level consortium qualification (expected during development)."
fi

echo
echo "Done. GDC references are installed in: $DEST"
echo "Normal Snakemake runs will use this installation and will not access the network."
