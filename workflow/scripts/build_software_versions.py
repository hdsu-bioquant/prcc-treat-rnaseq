#!/usr/bin/env python3
"""Build portable software-version provenance and MultiQC version metadata.

Production runs intentionally do not execute one version-probe job per tool.
Instead, the pipeline release owns ``workflow/config/software_versions.yaml`` as
its declared software/container manifest. This helper records the subset of
that manifest used by the current run and adds the actual Snakemake controller
version, which is inherently site/runtime-specific.

The same call writes:
  * results/run/software_versions.tsv
  * intermediate/qc/prcc_pipeline_mqc_versions.yml

The MultiQC YAML uses the supported ``*_mqc_versions.yml`` convention so that
MultiQC's native Software Versions section is complete and deterministic.
"""

from __future__ import annotations

import csv
from importlib import metadata as importlib_metadata
import json
import os
from pathlib import Path
from typing import Mapping, Sequence

import yaml


REQUIRED_TOOL_FIELDS = {"software", "version", "uri", "pull_default"}


def load_manifest(path: str | os.PathLike[str]) -> dict:
    with open(path) as fh:
        manifest = yaml.safe_load(fh)
    if not isinstance(manifest, dict) or not isinstance(manifest.get("tools"), dict):
        raise ValueError("software manifest must contain a top-level 'tools' mapping")
    for tool_id, spec in manifest["tools"].items():
        if not isinstance(spec, dict):
            raise ValueError(f"software manifest entry {tool_id!r} must be a mapping")
        missing = REQUIRED_TOOL_FIELDS - set(spec)
        if missing:
            raise ValueError(
                f"software manifest entry {tool_id!r} is missing: {', '.join(sorted(missing))}"
            )
    return manifest


def snakemake_version() -> str:
    try:
        return importlib_metadata.version("snakemake")
    except importlib_metadata.PackageNotFoundError:
        return "unknown"


def _resolved_display(value: str) -> str:
    if os.path.isabs(value):
        return os.path.basename(value)
    return value


def build_software_versions(
    manifest_path: str | os.PathLike[str],
    used_tools: Sequence[str],
    resolved_containers: Mapping[str, str],
    output_tsv: str | os.PathLike[str],
    output_multiqc_yml: str | os.PathLike[str],
) -> None:
    manifest = load_manifest(manifest_path)
    tools = manifest["tools"]

    unknown = sorted(set(used_tools) - set(tools))
    if unknown:
        raise ValueError("used tool(s) missing from software manifest: " + ", ".join(unknown))

    output_tsv = Path(output_tsv)
    output_multiqc_yml = Path(output_multiqc_yml)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    output_multiqc_yml.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for tool_id in sorted(used_tools):
        spec = tools[tool_id]
        resolved = resolved_containers.get(tool_id, spec["uri"])
        rows.append(
            {
                "tool_id": tool_id,
                "software": str(spec["software"]),
                "version": str(spec["version"]),
                "version_source": "pipeline_manifest",
                "container_source": str(spec["uri"]),
                "resolved_container": _resolved_display(str(resolved)),
            }
        )

    # Snakemake runs the workflow controller outside the per-rule containers, so
    # record the actual environment version instead of declaring it in the image manifest.
    rows.append(
        {
            "tool_id": "snakemake",
            "software": "Snakemake",
            "version": snakemake_version(),
            "version_source": "runtime",
            "container_source": "host/controller environment",
            "resolved_container": "host/controller environment",
        }
    )

    with output_tsv.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "tool_id",
                "software",
                "version",
                "version_source",
                "container_source",
                "resolved_container",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    # Flat software: version pairs produce the normal two-column MultiQC Software
    # Versions section. JSON-quoted scalars are valid YAML and prevent versions
    # such as 1.10 from being coerced to numbers.
    seen = set()
    with output_multiqc_yml.open("w") as fh:
        fh.write("# Generated from results/run/software_versions.tsv\n")
        for row in rows:
            software = row["software"]
            version = row["version"]
            if software in seen:
                raise ValueError(f"duplicate software display name: {software}")
            seen.add(software)
            fh.write(f"{json.dumps(software)}: {json.dumps(version)}\n")
