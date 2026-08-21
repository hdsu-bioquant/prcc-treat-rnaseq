#!/usr/bin/env python3
"""Build the canonical per-library gene-expression table.

The primary count is STAR's unstranded GeneCounts column for every assay. The two
stranded STAR columns are retained as diagnostics. GDC-style FPKM, FPKM-UQ and TPM
are computed from unstranded counts for full-length libraries only; QuantSeq gets
literal NA values because gene-length normalization is not appropriate for 3' tags.
UMI-bearing libraries additionally receive molecule-level HTSeq counts.

Usage:
  build_gene_expression.py <ReadsPerGene.out.tab> <gene_lengths.tsv> <assay> \
      <umi_counts.tsv|-> <out.tsv>
"""

import sys
import numpy as np
import pandas as pd

counts_f, lengths_f, assay, umi_f, out = sys.argv[1:6]
if assay not in {"full_length", "quantseq"}:
    raise SystemExit(f"Unsupported assay: {assay}")

rpg = pd.read_csv(
    counts_f,
    sep="\t",
    header=None,
    names=["gene_id", "unstranded", "stranded_first", "stranded_second"],
)
rpg["gene_id"] = rpg["gene_id"].astype(str)
rpg = rpg[~rpg["gene_id"].str.startswith("N_")].copy()
for col in ("unstranded", "stranded_first", "stranded_second"):
    rpg[col] = pd.to_numeric(rpg[col], errors="raise").astype(int)

lengths = pd.read_csv(lengths_f, sep="\t", dtype={"gene_id": str})
required = {"gene_id", "gene_length", "gene_type", "gene_name", "chromosome"}
missing = required.difference(lengths.columns)
if missing:
    raise SystemExit("gene_lengths.tsv is missing columns: " + ", ".join(sorted(missing)))

df = rpg.merge(lengths, on="gene_id", how="left", validate="one_to_one")
df["gene_length"] = pd.to_numeric(df["gene_length"], errors="coerce").fillna(0).astype(int)
df["gene_type"] = df["gene_type"].fillna("NA")
df["gene_name"] = df["gene_name"].fillna("NA")

if assay == "full_length":
    count = df["unstranded"].astype(float)
    length = df["gene_length"].replace(0, np.nan).astype(float)
    protein_coding = df["gene_type"] == "protein_coding"

    # Match the current NCI GDC gdc-rnaseq-tool definitions.  FPKM uses the
    # total unstranded count assigned to protein-coding genes.  FPKM-UQ uses
    # U = the 75th percentile of positive autosomal protein-coding counts and
    # G = the number of annotated autosomal protein-coding genes:
    #     FPKM-UQ = C * 1e9 / (U * G * L)
    rc_pc = float(count[protein_coding].sum())

    chromosome = df["chromosome"].fillna("NA").astype(str)
    non_autosomal = chromosome.isin({"chrX", "chrY", "chrM", "X", "Y", "M", "MT"})
    autosomal_pc = protein_coding & ~non_autosomal
    positive_autosomal_pc = autosomal_pc & (count > 0)
    uq_counts = count[positive_autosomal_pc]
    uq = float(np.quantile(uq_counts, 0.75)) if len(uq_counts) else 0.0
    n_autosomal_pc = int(autosomal_pc.sum())

    df["fpkm"] = (count * 1e9 / (rc_pc * length)).fillna(0.0) if rc_pc > 0 else 0.0
    if uq > 0 and n_autosomal_pc > 0:
        df["fpkm_uq"] = (count * 1e9 / (uq * n_autosomal_pc * length)).fillna(0.0)
    else:
        df["fpkm_uq"] = 0.0

    rate = count / length
    sum_rate = float(rate.sum(skipna=True))
    df["tpm"] = (rate / sum_rate * 1e6).fillna(0.0) if sum_rate > 0 else 0.0
else:
    df["fpkm"] = pd.NA
    df["fpkm_uq"] = pd.NA
    df["tpm"] = pd.NA

if umi_f != "-":
    umi = pd.read_csv(umi_f, sep="\t", header=None, names=["gene_id", "umi_molecule_count"])
    umi["gene_id"] = umi["gene_id"].astype(str)
    umi = umi[~umi["gene_id"].str.startswith("__")].copy()
    umi["umi_molecule_count"] = pd.to_numeric(umi["umi_molecule_count"], errors="raise").astype(int)
    df = df.merge(umi, on="gene_id", how="left", validate="one_to_one")
    df["umi_molecule_count"] = df["umi_molecule_count"].fillna(0).astype("Int64")
else:
    df["umi_molecule_count"] = pd.Series([pd.NA] * len(df), dtype="Int64")

columns = [
    "gene_id",
    "gene_name",
    "gene_type",
    "gene_length",
    "unstranded",
    "stranded_first",
    "stranded_second",
    "fpkm",
    "fpkm_uq",
    "tpm",
    "umi_molecule_count",
]
df[columns].to_csv(out, sep="\t", index=False, na_rep="NA", float_format="%.4f")
