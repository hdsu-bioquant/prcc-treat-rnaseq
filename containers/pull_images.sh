#!/usr/bin/env bash
#==============================================================================#
# Pre-pull the pipeline's container images to local .sif files, with retries,
# so the pipeline run never contacts a registry (avoids quay.io TLS timeouts).
# Container URIs are read from workflow/config/software_versions.yaml, the same
# manifest used by the workflow and run provenance.
#
# Run once on a node WITH internet (login node is fine), after activating the
# Snakemake environment (PyYAML is available there):
#   module load system/singularity
#   conda activate snakemake-9.19.0
#   bash containers/pull_images.sh
# Tunables: RETRIES (default 25), SLEEP seconds between tries (default 15),
#           ALL=1 to also pull the heavy optional-module images (fusion/TE/ASE).
#==============================================================================#
set -u
cd "$(dirname "$0")/.."          # -> RNA-Seq-pRCC
mkdir -p containers/sif
RETRIES="${RETRIES:-25}"; SLEEP="${SLEEP:-15}"
MANIFEST="workflow/config/software_versions.yaml"

rows_file="$(mktemp)"
trap 'rm -f "$rows_file"' EXIT

if ! python - "$MANIFEST" "${ALL:-0}" > "$rows_file" <<'PY'
import sys
try:
    import yaml
except ImportError as exc:
    raise SystemExit(
        "PyYAML is required to read workflow/config/software_versions.yaml; "
        "activate the Snakemake environment before running pull_images.sh"
    ) from exc

manifest_path, all_flag = sys.argv[1:3]
with open(manifest_path) as fh:
    manifest = yaml.safe_load(fh)
tools = manifest.get("tools", {}) if isinstance(manifest, dict) else {}
if not tools:
    raise SystemExit(f"No tools found in {manifest_path}")

pull_all = all_flag == "1"
for name, spec in tools.items():
    if pull_all or bool(spec.get("pull_default", False)):
        uri = spec.get("uri")
        if not uri:
            raise SystemExit(f"Missing uri for tool {name!r} in {manifest_path}")
        print(f"{name}\t{uri}")
PY
then
  echo "Failed to read container manifest: $MANIFEST" >&2
  exit 1
fi

fail=0
while IFS=$'\t' read -r name uri; do
  [[ -z "$name" ]] && continue
  out="containers/sif/${name}.sif"
  if [[ -s "$out" ]]; then echo "[skip] $name (already present)"; continue; fi
  ok=0
  for try in $(seq 1 "$RETRIES"); do
    echo "[$name] attempt $try/$RETRIES -> $uri"
    if apptainer pull --force "$out" "$uri"; then ok=1; echo "  [ok] $out"; break; fi
    echo "  ... failed; retry in ${SLEEP}s"; sleep "$SLEEP"
  done
  [[ $ok -ne 1 ]] && { echo "[FAIL] $name after $RETRIES tries"; fail=1; rm -f "$out"; }
done < "$rows_file"

echo "---"; ls -la containers/sif/ 2>/dev/null
[[ $fail -eq 0 ]] && echo "All requested images pulled." || { echo "Some images failed — re-run to resume (successful ones are skipped)."; exit 1; }
