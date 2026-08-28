#!/usr/bin/env bash
#==============================================================================#
# Pre-pull the pipeline's container images to local .sif files, with retries,
# so the pipeline run never contacts a registry (avoids quay.io TLS timeouts).
# Container URIs, optional local filenames, and optional version probes are read
# from workflow/config/software_versions.yaml, the same manifest used by the
# workflow and run provenance.
#
# Run once on a node WITH internet (login node is fine), after activating the
# Snakemake environment (PyYAML is available there):
#   module load system/singularity
#   conda activate prcc-rnaseq-controller
#   bash containers/pull_images.sh
# Tunables: RETRIES (default 25), SLEEP seconds between tries (default 15),
#           ALL=1 to also pull the heavy optional-module images (fusion/TE/ASE).
#==============================================================================#
set -u
cd "$(dirname "$0")/.."          # -> repository root
mkdir -p containers/sif
RETRIES="${RETRIES:-25}"; SLEEP="${SLEEP:-15}"
MANIFEST="workflow/config/software_versions.yaml"

if command -v apptainer >/dev/null 2>&1; then
  CONTAINER_RUNTIME="apptainer"
elif command -v singularity >/dev/null 2>&1; then
  CONTAINER_RUNTIME="singularity"
else
  echo "ERROR: neither apptainer nor singularity is available on PATH" >&2
  exit 1
fi

rows_file="$(mktemp)"
trap 'rm -f "$rows_file"' EXIT

if ! python - "$MANIFEST" "${ALL:-0}" > "$rows_file" <<'PY'
import os
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
        version = str(spec.get("version", ""))
        if not uri:
            raise SystemExit(f"Missing uri for tool {name!r} in {manifest_path}")
        local_sif = str(spec.get("local_sif", f"{name}.sif"))
        if os.path.basename(local_sif) != local_sif or local_sif in {"", ".", ".."}:
            raise SystemExit(
                f"Invalid local_sif for tool {name!r}: {local_sif!r}; expected a filename"
            )
        probe = str(spec.get("version_probe", ""))
        # Manifest values are maintainer-owned and may be interpreted by bash
        # below. Tabs/newlines would corrupt this compact transport format.
        for label, value in (("name", name), ("uri", uri), ("local_sif", local_sif),
                             ("version", version), ("version_probe", probe)):
            if "\t" in str(value) or "\n" in str(value):
                raise SystemExit(f"Unsupported tab/newline in {label} for tool {name!r}")
        print(f"{name}\t{uri}\t{local_sif}\t{version}\t{probe}")
PY
then
  echo "Failed to read container manifest: $MANIFEST" >&2
  exit 1
fi

verify_version() {
  local name="$1" image="$2" expected="$3" probe="$4"
  [[ -z "$probe" ]] && return 0

  local observed
  if ! observed="$("$CONTAINER_RUNTIME" exec "$image" bash -lc "$probe" 2>&1)"; then
    echo "  [verify-fail] $name: version probe failed: $probe" >&2
    echo "$observed" >&2
    return 1
  fi
  if [[ "$observed" != *"$expected"* ]]; then
    echo "  [verify-fail] $name: image did not report expected version $expected" >&2
    echo "  observed: $observed" >&2
    return 1
  fi
  echo "  [verified] $name reports version $expected"
  return 0
}

fail=0
while IFS=$'\t' read -r name uri local_sif version probe; do
  [[ -z "$name" ]] && continue
  out="containers/sif/${local_sif}"

  # The first 1.1.6 development upgrade used a temporary versioned filename to
  # prevent an old 1.1.4 umitools.sif from shadowing the new URI.  Now that
  # existing images are version-probed before reuse, migrate that verified image
  # back to the stable historical name and remove the temporary filename.
  if [[ "$name" == "umitools" && "$local_sif" == "umitools.sif" ]]; then
    migrated="containers/sif/umitools-1.1.6.sif"
    if [[ -s "$migrated" ]]; then
      if verify_version "$name" "$migrated" "$version" "$probe"; then
        if [[ -s "$out" ]] && verify_version "$name" "$out" "$version" "$probe"; then
          rm -f "$migrated"
          echo "[cleanup] removed redundant $migrated"
        else
          rm -f "$out"
          mv -f "$migrated" "$out"
          echo "[migrate] verified UMI-tools $version image -> $out"
        fi
      else
        echo "[cleanup] removing invalid temporary UMI-tools image: $migrated" >&2
        rm -f "$migrated"
      fi
    fi
  fi

  if [[ -s "$out" ]]; then
    if verify_version "$name" "$out" "$version" "$probe"; then
      echo "[skip] $name (already present: $out)"
      continue
    fi
    echo "[refresh] $name existing image does not match the maintained manifest"
  fi

  ok=0
  for try in $(seq 1 "$RETRIES"); do
    tmp="${out}.pulling.$$"
    rm -f "$tmp"
    echo "[$name] attempt $try/$RETRIES -> $uri"
    if "$CONTAINER_RUNTIME" pull --force "$tmp" "$uri"; then
      if verify_version "$name" "$tmp" "$version" "$probe"; then
        mv -f "$tmp" "$out"
        ok=1
        echo "  [ok] $out"
        break
      fi
    fi
    rm -f "$tmp"
    echo "  ... failed; retry in ${SLEEP}s"
    sleep "$SLEEP"
  done
  [[ $ok -ne 1 ]] && { echo "[FAIL] $name after $RETRIES tries"; fail=1; }
done < "$rows_file"

echo "---"; ls -la containers/sif/ 2>/dev/null
[[ $fail -eq 0 ]] && echo "All requested images pulled and verified." || { echo "Some images failed — re-run to resume (successful ones are skipped)."; exit 1; }
