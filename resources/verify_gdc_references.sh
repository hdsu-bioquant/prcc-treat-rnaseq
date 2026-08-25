#!/usr/bin/env bash
#==============================================================================#
# Verify an installed pRCC-RNA-Seq GDC reference bundle.
#
# Default verification is intentionally fast: it checks the expected FASTA/GTF
# and the essential files in the pre-built GDC STAR index. Normal workflow
# startup performs equivalent structural checks and does NOT hash ~30 GB of
# reference files on every run.
#
# Integrity layers:
#   --archives                  verify retained download archives against the
#                               official GDC MD5 values in gdc_resources.tsv
#   --canonical                verify the extracted installation against the
#                               maintainer-frozen canonical SHA256 manifest
#   --qualify                  perform --canonical and write/update the tiny
#                               site qualification stamp used by consortium runs
#   --write-canonical-manifest MAINTAINER ONLY: generate/replace the canonical
#                               installed-reference SHA256 manifest in resources/
#==============================================================================#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESOURCE_TABLE="$SCRIPT_DIR/gdc_resources.tsv"
CANONICAL_MANIFEST="$SCRIPT_DIR/gdc_installed_reference.sha256"
DEFAULT_DEST="$SCRIPT_DIR/gdc"
QUALIFICATION_STAMP=".prcc_treat_reference_qualification.tsv"

CHECK_ARCHIVES=false
CHECK_CANONICAL=false
QUALIFY=false
WRITE_CANONICAL_MANIFEST=false
DEST=""

usage() {
    cat <<EOF_USAGE
Usage: $0 [options] [REFERENCE_DIR]

Default behavior performs a fast structural check only.

Options:
  --archives                  Verify retained archives with official GDC MD5s.
  --canonical                Verify installed files against the maintainer-owned
                             resources/gdc_installed_reference.sha256 manifest.
  --qualify                  Verify against the canonical manifest and write the
                             site qualification stamp used by consortium runs.
  --write-canonical-manifest MAINTAINER ONLY: hash the qualified extracted bundle
                             and write resources/gdc_installed_reference.sha256.
  -h, --help                 Show this help.

Default REFERENCE_DIR: $DEFAULT_DEST
EOF_USAGE
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --archives)
            CHECK_ARCHIVES=true
            ;;
        --canonical)
            CHECK_CANONICAL=true
            ;;
        --qualify)
            CHECK_CANONICAL=true
            QUALIFY=true
            ;;
        --write-canonical-manifest)
            WRITE_CANONICAL_MANIFEST=true
            ;;
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
if [[ ! -d "$DEST" ]]; then
    echo "ERROR: reference directory does not exist: $DEST" >&2
    exit 1
fi
DEST="$(cd "$DEST" && pwd)"

if [[ ! -r "$RESOURCE_TABLE" ]]; then
    echo "ERROR: resource manifest not readable: $RESOURCE_TABLE" >&2
    exit 1
fi

check_file() {
    local path="$1"
    local label="$2"
    if [[ ! -s "$path" ]]; then
        echo "ERROR: missing or empty $label: $path" >&2
        return 1
    fi
}

check_star_index() {
    local index_dir="$1"
    if [[ ! -d "$index_dir" ]]; then
        echo "ERROR: missing STAR index directory: $index_dir" >&2
        return 1
    fi

    local required=(
        Genome
        SA
        SAindex
        chrLength.txt
        chrName.txt
        chrNameLength.txt
        chrStart.txt
        genomeParameters.txt
    )
    local f
    for f in "${required[@]}"; do
        check_file "$index_dir/$f" "STAR index file" || return 1
    done
}

GENOME_REL=""
GTF_REL=""
STAR_REL=""

while IFS=$'\t' read -r resource_id role url archive official_md5 installed_path; do
    [[ "$resource_id" == "resource_id" ]] && continue
    [[ -z "$resource_id" ]] && continue
    case "$role" in
        genome_fasta)
            GENOME_REL="$installed_path"
            check_file "$DEST/$installed_path" "genome FASTA"
            ;;
        annotation_gtf)
            GTF_REL="$installed_path"
            check_file "$DEST/$installed_path" "annotation GTF"
            ;;
        star_index)
            STAR_REL="$installed_path"
            check_star_index "$DEST/$installed_path"
            ;;
        *)
            echo "ERROR: unsupported role '$role' in $RESOURCE_TABLE" >&2
            exit 1
            ;;
    esac

    if [[ "$CHECK_ARCHIVES" == true ]]; then
        if [[ ! -s "$DEST/$archive" ]]; then
            echo "ERROR: archive missing for MD5 verification: $DEST/$archive" >&2
            exit 1
        fi
        echo "$official_md5  $DEST/$archive" | md5sum -c -
    fi
done < "$RESOURCE_TABLE"

if [[ -z "$GENOME_REL" || -z "$GTF_REL" || -z "$STAR_REL" ]]; then
    echo "ERROR: $RESOURCE_TABLE must define genome_fasta, annotation_gtf and star_index rows" >&2
    exit 1
fi

write_canonical_manifest() {
    local target="$CANONICAL_MANIFEST"
    local tmp="${target}.tmp.$$"
    rm -f "$tmp"

    (
        cd "$DEST"
        {
            printf '%s\0' "$GENOME_REL" "$GTF_REL"
            find "$STAR_REL" -type f -print0
        } | sort -z | while IFS= read -r -d '' path; do
            sha256sum "$path"
        done
    ) > "$tmp"

    mv "$tmp" "$target"
    echo "Wrote canonical installed-reference manifest: $target"
    echo "IMPORTANT: this file is maintainer-owned. Freeze/commit it only after the"
    echo "           reference installation has passed the realistic qualification."
}

verify_canonical_manifest() {
    if [[ ! -s "$CANONICAL_MANIFEST" ]]; then
        echo "ERROR: canonical installed-reference manifest is not available:" >&2
        echo "  $CANONICAL_MANIFEST" >&2
        echo "It should be generated and frozen by maintainers only after realistic" >&2
        echo "reference qualification; partners should not create their own baseline." >&2
        exit 1
    fi

    (
        cd "$DEST"
        sha256sum --check "$CANONICAL_MANIFEST"
    )
}

write_qualification_stamp() {
    local manifest_sha256
    manifest_sha256="$(sha256sum "$CANONICAL_MANIFEST" | awk '{print $1}')"
    local target="$DEST/$QUALIFICATION_STAMP"
    local tmp="${target}.tmp.$$"

    {
        printf 'canonical_manifest_sha256\t%s\n' "$manifest_sha256"
        printf 'qualified_at_utc\t%s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')"
    } > "$tmp"
    mv "$tmp" "$target"
    echo "Wrote consortium qualification stamp: $target"
}

if [[ "$WRITE_CANONICAL_MANIFEST" == true ]]; then
    write_canonical_manifest
fi

if [[ "$CHECK_CANONICAL" == true ]]; then
    verify_canonical_manifest
fi

if [[ "$QUALIFY" == true ]]; then
    write_qualification_stamp
fi

echo "PASS: GDC reference installation is structurally complete: $DEST"
