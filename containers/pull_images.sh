#!/usr/bin/env bash
#==============================================================================#
# Pre-pull the pipeline's container images to local .sif files, with retries,
# so the pipeline run never contacts a registry (avoids quay.io TLS timeouts).
# common.smk automatically uses containers/sif/<name>.sif when present.
#
# Run once on a node WITH internet (login node is fine):
#   module load system/singularity
#   bash containers/pull_images.sh
# Tunables: RETRIES (default 25), SLEEP seconds between tries (default 15),
#           ALL=1 to also pull the heavy optional-module images (fusion/TE/ASE/DE).
#==============================================================================#
set -u
cd "$(dirname "$0")/.."          # -> RNA-Seq-pRCC
mkdir -p containers/sif
RETRIES="${RETRIES:-25}"; SLEEP="${SLEEP:-15}"

# core images used by the validation / default pipeline (full-length + QuantSeq)
names=(py star fastqc multiqc fastp bbmap samtools htseq umitools)
uris=(
  "docker://quay.io/biocontainers/pandas:1.5.2"
  "docker://quay.io/biocontainers/star:2.7.5c--0"
  "docker://quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
  "docker://quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0"
  "docker://quay.io/biocontainers/fastp:0.23.4--hadf994f_2"
  "docker://quay.io/biocontainers/bbmap:39.06--h92535d8_0"
  "docker://quay.io/biocontainers/samtools:1.19--h50ea8bc_0"
  "docker://quay.io/biocontainers/htseq:2.0.9--py39h918f1d6_0"
  "docker://quay.io/biocontainers/umi_tools:1.1.4--py39hf95cd2a_2"
)

if [[ "${ALL:-0}" == "1" ]]; then
  names+=(rseqc arriba starfusion tetx gatk de)
  uris+=(
    "docker://quay.io/biocontainers/rseqc:5.0.3--py39hf95cd2a_0"
    "docker://quay.io/biocontainers/arriba:1.1.0--h2e03b76_3"
    "docker://trinityctat/starfusion:1.6.0"
    "docker://quay.io/biocontainers/tetranscripts:2.2.3--pyhdfd78af_0"
    "docker://broadinstitute/gatk:4.5.0.0"
    "docker://bioconductor/bioconductor_docker:RELEASE_3_18"
  )
fi

fail=0
for i in "${!names[@]}"; do
  name="${names[$i]}"; uri="${uris[$i]}"; out="containers/sif/${name}.sif"
  if [[ -s "$out" ]]; then echo "[skip] $name (already present)"; continue; fi
  ok=0
  for try in $(seq 1 "$RETRIES"); do
    echo "[$name] attempt $try/$RETRIES -> $uri"
    if apptainer pull --force "$out" "$uri"; then ok=1; echo "  [ok] $out"; break; fi
    echo "  ... failed; retry in ${SLEEP}s"; sleep "$SLEEP"
  done
  [[ $ok -ne 1 ]] && { echo "[FAIL] $name after $RETRIES tries"; fail=1; rm -f "$out"; }
done

echo "---"; ls -la containers/sif/ 2>/dev/null
[[ $fail -eq 0 ]] && echo "All requested images pulled." || { echo "Some images failed — re-run to resume (successful ones are skipped)."; exit 1; }
