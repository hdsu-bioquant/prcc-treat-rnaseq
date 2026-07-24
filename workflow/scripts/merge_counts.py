#!/usr/bin/env python
"""Cohort raw-count matrix (gene x sample) from per-sample STAR gene-count tables.

Uses the UNSTRANDED STAR column (the pRCC-TREAT primary) for BOTH branches, so the
matrix is uniform (full-length + QuantSeq on the same non-deduplicated STAR basis).
`assay` is carried as a header annotation row for downstream covariate modeling.
Usage: merge_counts.py <samples.tsv> <results_dir> <out.tsv>
"""
import sys, os
import pandas as pd

samples_f, results, out = sys.argv[1], sys.argv[2], sys.argv[3]
samples = pd.read_csv(samples_f, sep="\t", dtype=str).fillna("-")

series, assays = {}, {}
for _, r in samples.iterrows():
    s, assay = r["sample"], r["assay"]
    branch = "full_length" if assay == "full_length_pe" else "quantseq"
    f = os.path.join(results, branch, s, s + ".star_gene_counts.tsv")
    d = pd.read_csv(f, sep="\t")
    series[s] = d.set_index("gene_id")["unstranded"]
    assays[s] = assay

mat = pd.DataFrame(series).fillna(0).astype(int)
mat.index.name = "gene_id"
with open(out, "w") as o:
    o.write("# assay\t" + "\t".join(assays[c] for c in mat.columns) + "\n")
mat.to_csv(out, sep="\t", mode="a")
