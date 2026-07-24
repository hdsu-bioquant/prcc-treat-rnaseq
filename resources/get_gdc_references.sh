#!/usr/bin/env bash
#==============================================================================#
# Download the EXACT GDC reference resources (GRCh38.d1.vd1 + GENCODE v36 +
# the GDC-built STAR 2.7.5c index) and verify MD5s. This is the alternative to
# letting the Snakemake references.smk rules fetch them.
#
# Source: https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files
#==============================================================================#
set -euo pipefail
DEST="${1:-resources/gdc}"
mkdir -p "$DEST"; cd "$DEST"

dl () { # url  outfile  md5
  echo ">> $2"; wget -O "$2" "$1"
  echo "$3  $2" | md5sum -c -
}

# Genome FASTA (GRCh38.d1.vd1)
dl "https://api.gdc.cancer.gov/data/254f697d-310d-4d7d-a27b-27fbf767a834" \
   "GRCh38.d1.vd1.fa.tar.gz" "3ffbcfe2d05d43206f57f81ebb251dc9"
tar -xzf GRCh38.d1.vd1.fa.tar.gz

# GENCODE v36 GTF
dl "https://api.gdc.cancer.gov/data/be002a2c-3b27-43f3-9e0f-fd47db92a6b5" \
   "gencode.v36.annotation.gtf.gz" "c03931958d4572148650d62eb6dec41a"
gunzip -f gencode.v36.annotation.gtf.gz

# STAR 2.7.5c index (sjdbOverhang 100), ~25 GB
dl "https://api.gdc.cancer.gov/data/c0008693-0583-4eac-bd5c-583070763893" \
   "star-2.7.5c_GRCh38.d1.vd1_gencode.v36.tgz" "acafb76bba5e3e80eb028dc05f002ffc"
mkdir -p star-2.7.5c_GRCh38.d1.vd1_gencode.v36
tar -xzf star-2.7.5c_GRCh38.d1.vd1_gencode.v36.tgz \
    -C star-2.7.5c_GRCh38.d1.vd1_gencode.v36 --strip-components=1

echo "Done. References in $DEST"
