#!/usr/bin/env python
"""Augment STAR ReadsPerGene.out.tab with FPKM, FPKM-UQ and TPM, GDC-style.

Reproduces the GDC "rna_seq.augmented_star_gene_counts" output. STAR's
ReadsPerGene.out.tab has 4 header rows (N_unmapped, N_multimapping, N_noFeature,
N_ambiguous) then per gene: gene_id, unstranded, stranded_first, stranded_second.
GDC uses the UNSTRANDED column as the headline and computes:
  FPKM    = RCg * 1e9 / (RCpc * L)          RCpc = total reads over protein-coding genes
  FPKM-UQ = RCg * 1e9 / (RC75 * L)          RC75 = 75th-percentile protein-coding gene count
  TPM     = (RCg/L) / sum_g(RCg/L) * 1e6
Ref: https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/
     https://docs.gdc.cancer.gov/Encyclopedia/pages/FPKM-UQ/  (FPKM, FPKM-UQ)
NOTE: validate against a GDC reference sample before production use.
Usage: augment_star_counts.py <ReadsPerGene.out.tab> <gene_lengths.tsv> <col> <out.tsv>
  <col> = unstranded | forward | reverse
"""
import sys
import pandas as pd
import numpy as np

counts_f, lengths_f, col, out = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
colmap = {"unstranded": 1, "forward": 2, "reverse": 3}
ci = colmap[col]

rpg = pd.read_csv(counts_f, sep="\t", header=None,
                  names=["gene_id", "unstranded", "stranded_first", "stranded_second"])
rpg = rpg[~rpg["gene_id"].str.startswith("N_")].copy()        # drop the 4 summary rows

ln = pd.read_csv(lengths_f, sep="\t")
df = rpg.merge(ln, on="gene_id", how="left")
df["gene_length"] = df["gene_length"].fillna(0).astype(float)

count = df.iloc[:, ci].astype(float)                          # chosen strand column
df["count"] = count

pc = df["gene_type"] == "protein_coding"
RCpc = float(count[pc].sum())
pc_counts = count[pc]
RC75 = float(np.percentile(pc_counts, 75)) if len(pc_counts) else 0.0

L = df["gene_length"].replace(0, np.nan)
df["fpkm_unstranded"]    = (count * 1e9 / (RCpc * L)).fillna(0) if RCpc > 0 else 0.0
df["fpkm_uq_unstranded"] = (count * 1e9 / (RC75 * L)).fillna(0) if RC75 > 0 else 0.0
rate = (count / L)
sum_rate = float(rate.sum(skipna=True))
df["tpm_unstranded"] = (rate / sum_rate * 1e6).fillna(0) if sum_rate > 0 else 0.0

cols = ["gene_id", "unstranded", "stranded_first", "stranded_second",
        "gene_length", "gene_type", "gene_name",
        "fpkm_unstranded", "fpkm_uq_unstranded", "tpm_unstranded", "count"]
df[cols].to_csv(out, sep="\t", index=False)
