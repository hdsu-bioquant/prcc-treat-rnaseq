#!/usr/bin/env python3
"""Convert canonical run QC metrics to a MultiQC 1.21-compatible custom table.

The canonical machine-readable QC contract is ``results/qc/qc_metrics.tsv``.
This helper creates an embedded-config ``*_mqc.tsv`` file with one row per
sequencing library.  Using TSV here is deliberate: it exercises MultiQC's
long-standing custom-table parser and avoids version-specific ambiguity in the
stand-alone JSON table parser used by MultiQC 1.21.
"""

import csv
import sys

import pandas as pd

inp, out = sys.argv[1:3]
df = pd.read_csv(inp, sep="\t", na_values=["NA"], keep_default_na=True)

required = [
    "library_id",
    "sample_id",
    "assay",
    "layout",
    "has_umi",
    "star_input_records",
    "uniquely_mapped_percent",
    "multi_mapped_percent",
    "gene_assigned_unstranded",
    "gene_assigned_percent",
    "umi_molecules_assigned",
]
missing = [c for c in required if c not in df.columns]
if missing:
    raise SystemExit("QC table is missing required column(s): " + ", ".join(missing))
if df["library_id"].isna().any() or df["library_id"].duplicated().any():
    raise SystemExit("QC table must contain one unique, non-empty library_id per row")


def assay_label(value):
    return {"full_length": "Full-length", "quantseq": "QuantSeq"}.get(str(value), str(value))


def layout_label(value):
    return {"paired": "Paired-end", "single": "Single-end"}.get(str(value), str(value))


def yes_no(value):
    return "Yes" if str(value).lower() == "true" else "No"


def text_value(value):
    if pd.isna(value):
        return ""
    if hasattr(value, "item"):
        value = value.item()
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value)


# MultiQC custom-content TSV files accept YAML configuration in leading comment
# lines.  For table plots, column configuration belongs in the top-level
# ``headers`` mapping (not inside ``pconfig``).
meta = """# id: prcc_rnaseq_qc
# section_name: Library QC Summary
# description: >-
#   One row per sequencing library. STAR metrics are library-level. Detailed
#   FastQC R1/R2 diagnostics remain available in the dedicated FastQC sections
#   below. No consortium PASS/WARN/FAIL thresholds are applied yet.
# format: tsv
# plot_type: table
# pconfig:
#   id: prcc_rnaseq_qc_table
#   title: pRCC-RNA-Seq library QC summary
#   namespace: pRCC-RNA-Seq
#   col1_header: Library ID
# headers:
#   sample_id:
#     title: Sample ID
#     description: Biological sample identifier associated with this sequencing library.
#   assay:
#     title: Assay
#     description: Library assay type.
#   layout:
#     title: Layout
#     description: Sequencing read layout.
#   has_umi:
#     title: UMI
#     description: Whether the library contains a UMI used for molecule deduplication.
#   star_input_records:
#     title: STAR input
#     description: Records entering STAR after preprocessing.
#     format: '{:,.0f}'
#   uniquely_mapped_percent:
#     title: Unique mapped %
#     description: Percentage of STAR input records mapped uniquely.
#     suffix: '%'
#     format: '{:.1f}'
#     min: 0
#     max: 100
#   multi_mapped_percent:
#     title: Multi-mapped %
#     description: Percentage of STAR input records mapped to multiple loci.
#     suffix: '%'
#     format: '{:.1f}'
#     min: 0
#     max: 100
#   gene_assigned_unstranded:
#     title: Gene assigned
#     description: Sum of canonical unstranded STAR gene counts.
#     format: '{:,.0f}'
#   gene_assigned_percent:
#     title: Assigned %
#     description: Canonical unstranded gene-assigned count divided by STAR input records.
#     suffix: '%'
#     format: '{:.1f}'
#     min: 0
#     max: 100
#   umi_molecules_assigned:
#     title: UMI molecules
#     description: Assigned gene-level molecules after UMI deduplication; blank for non-UMI libraries.
#     format: '{:,.0f}'
"""

columns = [
    "library_id",
    "sample_id",
    "assay",
    "layout",
    "has_umi",
    "star_input_records",
    "uniquely_mapped_percent",
    "multi_mapped_percent",
    "gene_assigned_unstranded",
    "gene_assigned_percent",
    "umi_molecules_assigned",
]

with open(out, "w", newline="") as fh:
    fh.write(meta)
    writer = csv.writer(fh, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_MINIMAL)
    writer.writerow(columns)
    for _, row in df.iterrows():
        writer.writerow(
            [
                text_value(row["library_id"]),
                text_value(row["sample_id"]),
                assay_label(row["assay"]),
                layout_label(row["layout"]),
                yes_no(row["has_umi"]),
                text_value(row["star_input_records"]),
                text_value(row["uniquely_mapped_percent"]),
                text_value(row["multi_mapped_percent"]),
                text_value(row["gene_assigned_unstranded"]),
                text_value(row["gene_assigned_percent"]),
                text_value(row["umi_molecules_assigned"]),
            ]
        )
