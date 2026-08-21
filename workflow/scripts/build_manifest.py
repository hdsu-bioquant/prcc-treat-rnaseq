#!/usr/bin/env python3
"""Create portable-results manifest and checksum files.

The package checksum covers every portable file except the checksum files themselves.
The validation checksum intentionally covers only deterministic canonical data products;
it excludes MultiQC HTML, provenance timestamps, and environment-dependent metadata.

Usage:
  build_manifest.py <results_dir> <manifest.tsv> <checksums.sha256> <validation.sha256>
"""

from __future__ import annotations

import hashlib
import os
import sys
from pathlib import Path

results = Path(sys.argv[1]).resolve()
manifest = Path(sys.argv[2]).resolve()
checksums = Path(sys.argv[3]).resolve()
validation = Path(sys.argv[4]).resolve()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def rel(path: Path) -> str:
    return path.relative_to(results).as_posix()


def is_validation_file(r: str) -> bool:
    if r.startswith("libraries/") and (r.endswith("/gene_expression.tsv") or r.endswith("/qc_metrics.tsv")):
        return True
    if r.startswith("matrices/") and r.endswith(".tsv"):
        return True
    if r == "qc/qc_metrics.tsv":
        return True
    if r in {"run/libraries.tsv", "run/config.yaml", "run/references.tsv"}:
        return True
    return False


def file_type(r: str) -> str:
    if r.endswith("/gene_expression.tsv"):
        return "gene_expression"
    if r.endswith("/qc_metrics.tsv") or r == "qc/qc_metrics.tsv":
        return "qc_metrics"
    if r.startswith("matrices/"):
        return "matrix"
    if r == "qc/multiqc_report.html":
        return "multiqc_report"
    if r.startswith("run/"):
        return "run_metadata"
    return "result"


for p in (manifest, checksums, validation):
    p.parent.mkdir(parents=True, exist_ok=True)
    if p.exists():
        p.unlink()

excluded = {manifest, checksums, validation}
files = sorted(p for p in results.rglob("*") if p.is_file() and p.resolve() not in excluded)

with manifest.open("w") as out:
    out.write("path\tfile_type\tsize_bytes\tsha256\n")
    for path in files:
        r = rel(path)
        out.write(f"{r}\t{file_type(r)}\t{path.stat().st_size}\t{sha256(path)}\n")

with validation.open("w") as out:
    for path in files:
        r = rel(path)
        if is_validation_file(r):
            out.write(f"{sha256(path)}  {r}\n")

# Package integrity includes both manifests, but not the package checksum itself.
package_files = sorted(files + [manifest, validation])
with checksums.open("w") as out:
    for path in package_files:
        out.write(f"{sha256(path)}  {rel(path)}\n")
