#!/usr/bin/env python3
"""Load and validate the pRCC-TREAT sequencing-library sample sheet.

One row represents one sequencing library. ``library_id`` is the workflow/output
identifier; ``sample_id`` identifies the biological sample and may therefore occur in
multiple rows (for example, when a biological sample has multiple libraries).

The workflow currently supports the library layouts used by the consortium:
  - full_length + paired
  - quantseq + single

UMI metadata is library-specific and independent of assay. ``has_umi`` is always
required. The three UMI-detail columns are required only when at least one UMI-bearing
library is present; an all-non-UMI sheet may omit them entirely.
"""

from __future__ import annotations

import csv
import os
import re
from pathlib import Path
from typing import Iterable

import pandas as pd

from umi import UMI_COLUMNS, parse_umi_spec

BASE_REQUIRED_COLUMNS = [
    "library_id",
    "sample_id",
    "assay",
    "layout",
    "strandedness",
    "fq1",
    "fq2",
    "has_umi",
]

ALLOWED_ASSAYS = {"full_length", "quantseq"}
ALLOWED_LAYOUTS = {"paired", "single"}
ALLOWED_STRANDEDNESS = {"unstranded", "forward", "reverse"}
SUPPORTED_ASSAY_LAYOUTS = {
    ("full_length", "paired"),
    ("quantseq", "single"),
}
FASTQ_SUFFIXES = (".fastq", ".fastq.gz", ".fq", ".fq.gz")
ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
MISSING = {"", "-"}


def _duplicate_items(values: Iterable[str]) -> list[str]:
    seen = set()
    duplicates = []
    for value in values:
        if value in seen and value not in duplicates:
            duplicates.append(value)
        seen.add(value)
    return duplicates


def _check_header(path: Path) -> list[str]:
    with path.open("r", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            raise ValueError(f"Sample sheet is empty: {path}")

    if not header:
        raise ValueError(f"Sample sheet has no header: {path}")

    duplicates = _duplicate_items(header)
    if duplicates:
        raise ValueError(
            "Sample-sheet validation failed:\n  - duplicate column name(s): "
            + ", ".join(duplicates)
        )
    return header


def load_and_validate_samples(path: str | os.PathLike[str]) -> pd.DataFrame:
    """Return a normalized, validated sample-sheet DataFrame indexed by library_id.

    Extra metadata columns are allowed and preserved. Required technical columns are
    normalized to the canonical lower-case vocabulary where appropriate. Missing UMI
    columns are added internally as ``-`` for all-non-UMI sheets so downstream workflow
    code always receives one normalized schema.
    """

    path = Path(path)
    if not path.is_file():
        raise ValueError(f"Sample sheet does not exist or is not a file: {path}")

    header = _check_header(path)
    errors: list[str] = []

    missing_base_columns = [c for c in BASE_REQUIRED_COLUMNS if c not in header]
    if missing_base_columns:
        errors.append("missing required column(s): " + ", ".join(missing_base_columns))

    unexpected_blank_headers = [i + 1 for i, name in enumerate(header) if not name.strip()]
    if unexpected_blank_headers:
        errors.append(
            "blank column name(s) at position(s): "
            + ", ".join(map(str, unexpected_blank_headers))
        )

    if errors:
        raise ValueError("Sample-sheet validation failed:\n  - " + "\n  - ".join(errors))

    df = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        na_filter=False,
    )

    if df.empty:
        raise ValueError("Sample-sheet validation failed:\n  - sample sheet contains no libraries")

    # Strip whitespace from all fields without altering the user's extra columns.
    for col in df.columns:
        df[col] = df[col].astype(str).str.strip()

    # Canonicalize workflow-driving categorical fields that are always present.
    for col in ("assay", "layout", "strandedness", "has_umi"):
        df[col] = df[col].str.lower()

    # UMI-detail columns are conditionally required at the input interface. Determine
    # whether they are needed before normalizing an all-non-UMI sheet to the full schema.
    has_any_umi = bool((df["has_umi"] == "true").any())
    present_umi_columns = [c for c in UMI_COLUMNS if c in header]
    missing_umi_columns = [c for c in UMI_COLUMNS if c not in header]
    if present_umi_columns and missing_umi_columns:
        errors.append(
            "UMI-detail columns must be supplied together; missing: "
            + ", ".join(missing_umi_columns)
        )
    elif has_any_umi and missing_umi_columns:
        errors.append(
            "UMI-bearing libraries are present but UMI-detail columns are missing: "
            + ", ".join(missing_umi_columns)
        )

    for col in UMI_COLUMNS:
        if col not in df.columns:
            df[col] = "-"

    df["umi_pattern"] = df["umi_pattern"].str.upper()
    df["umi_location"] = df["umi_location"].str.lower()
    df.loc[df["fq2"] == "", "fq2"] = "-"
    for col in UMI_COLUMNS:
        df.loc[df[col] == "", col] = "-"

    # Required identifiers.
    for rownum, row in df.iterrows():
        line = rownum + 2  # header is line 1
        library_id = row["library_id"]
        sample_id = row["sample_id"]

        if not library_id:
            errors.append(f"line {line}: library_id is empty")
        elif not ID_RE.fullmatch(library_id):
            errors.append(
                f"line {line}: library_id '{library_id}' contains unsupported characters; "
                "use letters, numbers, '.', '_' or '-' and do not start with punctuation"
            )

        if not sample_id:
            errors.append(f"line {line}: sample_id is empty")
        elif not ID_RE.fullmatch(sample_id):
            errors.append(
                f"line {line}: sample_id '{sample_id}' contains unsupported characters; "
                "use letters, numbers, '.', '_' or '-' and do not start with punctuation"
            )

    duplicate_libraries = df.loc[df["library_id"].duplicated(keep=False), "library_id"].unique()
    if len(duplicate_libraries):
        errors.append("library_id must be unique; duplicate(s): " + ", ".join(duplicate_libraries))

    # Controlled vocabulary, FASTQ existence and UMI semantics.
    for rownum, row in df.iterrows():
        line = rownum + 2
        library_id = row["library_id"] or f"line_{line}"
        assay = row["assay"]
        layout = row["layout"]
        strandedness = row["strandedness"]
        has_umi = row["has_umi"]

        if assay not in ALLOWED_ASSAYS:
            errors.append(
                f"line {line} ({library_id}): assay='{assay}' is unsupported; "
                f"expected one of {sorted(ALLOWED_ASSAYS)}"
            )
        if layout not in ALLOWED_LAYOUTS:
            errors.append(
                f"line {line} ({library_id}): layout='{layout}' is unsupported; "
                f"expected one of {sorted(ALLOWED_LAYOUTS)}"
            )
        if strandedness not in ALLOWED_STRANDEDNESS:
            errors.append(
                f"line {line} ({library_id}): strandedness='{strandedness}' is unsupported; "
                f"expected one of {sorted(ALLOWED_STRANDEDNESS)}"
            )
        if (
            assay in ALLOWED_ASSAYS
            and layout in ALLOWED_LAYOUTS
            and (assay, layout) not in SUPPORTED_ASSAY_LAYOUTS
        ):
            errors.append(
                f"line {line} ({library_id}): assay/layout combination '{assay} + {layout}' "
                "is not currently implemented (supported: full_length+paired, quantseq+single)"
            )
        if has_umi not in {"true", "false"}:
            errors.append(
                f"line {line} ({library_id}): has_umi='{has_umi}' is invalid; use true or false"
            )

        fq1 = row["fq1"]
        fq2 = row["fq2"]
        if fq1 in MISSING:
            errors.append(f"line {line} ({library_id}): fq1 is required")
        elif not fq1.endswith(FASTQ_SUFFIXES):
            errors.append(
                f"line {line} ({library_id}): fq1 '{fq1}' does not have a supported FASTQ suffix"
            )
        elif not Path(fq1).is_file():
            errors.append(f"line {line} ({library_id}): fq1 does not exist: {fq1}")

        if layout == "paired":
            if fq2 in MISSING:
                errors.append(f"line {line} ({library_id}): paired library requires fq2")
            elif not fq2.endswith(FASTQ_SUFFIXES):
                errors.append(
                    f"line {line} ({library_id}): fq2 '{fq2}' does not have a supported FASTQ suffix"
                )
            elif not Path(fq2).is_file():
                errors.append(f"line {line} ({library_id}): fq2 does not exist: {fq2}")
        elif layout == "single" and fq2 not in MISSING:
            errors.append(
                f"line {line} ({library_id}): single-end library must use '-' or an empty fq2"
            )

        umi_values = {col: row[col] for col in UMI_COLUMNS}
        if has_umi == "true":
            missing_values = [col for col, value in umi_values.items() if value in MISSING]
            if missing_values:
                errors.append(
                    f"line {line} ({library_id}): has_umi=true requires value(s) for "
                    + ", ".join(missing_values)
                )
            else:
                try:
                    spec = parse_umi_spec(
                        row["umi_pattern"], row["umi_location"], row["umi_discard_bases"]
                    )
                except ValueError as exc:
                    errors.append(f"line {line} ({library_id}): {exc}")
                else:
                    if spec.location == "read2_start" and layout != "paired":
                        errors.append(
                            f"line {line} ({library_id}): umi_location=read2_start requires layout=paired"
                        )
        elif has_umi == "false":
            supplied = {col: value for col, value in umi_values.items() if value not in MISSING}
            if supplied:
                rendered = ", ".join(f"{col}='{value}'" for col, value in supplied.items())
                errors.append(
                    f"line {line} ({library_id}): has_umi=false but UMI metadata was supplied "
                    f"({rendered}); use '-' or empty values"
                )

    # Detect accidental FASTQ reuse across libraries, including via symlinks.
    fastq_owners: dict[str, list[str]] = {}
    for _, row in df.iterrows():
        for col in ("fq1", "fq2"):
            path_value = row[col]
            if path_value in MISSING:
                continue
            canonical = os.path.realpath(os.path.abspath(path_value))
            fastq_owners.setdefault(canonical, []).append(f"{row['library_id']}:{col}")
    reused = {p: owners for p, owners in fastq_owners.items() if len(owners) > 1}
    for path_value, owners in reused.items():
        errors.append(
            "FASTQ is assigned more than once (possible sample-sheet mistake): "
            f"{path_value} -> {', '.join(owners)}"
        )

    if errors:
        raise ValueError("Sample-sheet validation failed:\n  - " + "\n  - ".join(errors))

    return df.set_index("library_id", drop=False)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Validate a pRCC-TREAT library sample sheet")
    parser.add_argument("samples_tsv")
    args = parser.parse_args()
    validated = load_and_validate_samples(args.samples_tsv)
    print(
        f"PASS: {len(validated)} libraries, "
        f"{validated['sample_id'].nunique()} biological samples"
    )
