#!/usr/bin/env python3
"""Build the stable one-row pRCC-RNA-Seq QC table for one library.

This deliberately exposes a small, version-controlled QC schema instead of making
MultiQC's internal export formats the consortium data contract. MultiQC remains the
human-facing report and consumes the run-level version of this table as custom content.

Usage:
  build_qc_metrics.py <star.Log.final.out> <gene_expression.tsv> <library_id> \
      <sample_id> <assay> <layout> <has_umi> <out.tsv>
"""

import sys
import pandas as pd

log_f, expr_f, library_id, sample_id, assay, layout, has_umi, out = sys.argv[1:9]


def parse_star_log(path):
    values = {}
    with open(path) as fh:
        for line in fh:
            if "|" not in line:
                continue
            key, value = line.split("|", 1)
            values[key.strip()] = value.strip()
    return values


def as_int(value):
    if value is None or value == "":
        return pd.NA
    return int(value.replace(",", ""))


def as_percent(value):
    if value is None or value == "":
        return pd.NA
    return float(value.rstrip("%"))


star = parse_star_log(log_f)
expr = pd.read_csv(expr_f, sep="\t", na_values=["NA"], keep_default_na=True)

star_input = as_int(star.get("Number of input reads"))
assigned = int(pd.to_numeric(expr["unstranded"], errors="raise").sum())
assigned_pct = (assigned / star_input * 100.0) if (not pd.isna(star_input) and star_input != 0) else pd.NA

if has_umi.lower() == "true":
    umi_col = pd.to_numeric(expr["umi_molecule_count"], errors="coerce")
    if umi_col.isna().any():
        raise SystemExit("UMI-bearing library has NA values in umi_molecule_count")
    umi_molecules = int(umi_col.sum())
else:
    umi_molecules = pd.NA

row = {
    "library_id": library_id,
    "sample_id": sample_id,
    "assay": assay,
    "layout": layout,
    "has_umi": has_umi.lower(),
    "star_input_records": star_input,
    "star_input_unit": "fragments" if layout == "paired" else "reads",
    "uniquely_mapped_reads": as_int(star.get("Uniquely mapped reads number")),
    "uniquely_mapped_percent": as_percent(star.get("Uniquely mapped reads %")),
    "multi_mapped_reads": as_int(star.get("Number of reads mapped to multiple loci")),
    "multi_mapped_percent": as_percent(star.get("% of reads mapped to multiple loci")),
    "too_many_loci_reads": as_int(star.get("Number of reads mapped to too many loci")),
    "too_many_loci_percent": as_percent(star.get("% of reads mapped to too many loci")),
    "unmapped_mismatches_percent": as_percent(star.get("% of reads unmapped: too many mismatches")),
    "unmapped_too_short_percent": as_percent(star.get("% of reads unmapped: too short")),
    "unmapped_other_percent": as_percent(star.get("% of reads unmapped: other")),
    "gene_assigned_unstranded": assigned,
    "gene_assigned_percent": assigned_pct,
    "umi_molecules_assigned": umi_molecules,
}

pd.DataFrame([row]).to_csv(out, sep="\t", index=False, na_rep="NA", float_format="%.6g")
