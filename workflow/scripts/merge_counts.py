#!/usr/bin/env python3
"""Cohort raw-count matrix (gene x library) from per-library STAR count tables.

Uses the UNSTRANDED STAR column (the pRCC-TREAT primary) for both assays. Matrix
columns are ``library_id`` values. Comment rows retain the biological ``sample_id``
and technical assay/layout mapping for downstream provenance.

Usage: merge_counts.py <samples.tsv> <results_dir> <out.tsv>
"""

import os
import sys

import pandas as pd


samples_f, results, out = sys.argv[1], sys.argv[2], sys.argv[3]
samples = pd.read_csv(samples_f, sep="\t", dtype=str, keep_default_na=False).set_index("library_id", drop=False)

series = {}
for library_id, row in samples.iterrows():
    assay = row["assay"]
    branch = "full_length" if assay == "full_length" else "quantseq"
    path = os.path.join(results, branch, library_id, library_id + ".star_gene_counts.tsv")
    counts = pd.read_csv(path, sep="\t")
    series[library_id] = counts.set_index("gene_id")["unstranded"]

matrix = pd.DataFrame(series).fillna(0).astype(int)
matrix.index.name = "gene_id"

with open(out, "w") as fh:
    ordered = list(matrix.columns)
    fh.write("# sample_id\t" + "\t".join(samples.loc[x, "sample_id"] for x in ordered) + "\n")
    fh.write("# assay\t" + "\t".join(samples.loc[x, "assay"] for x in ordered) + "\n")
    fh.write("# layout\t" + "\t".join(samples.loc[x, "layout"] for x in ordered) + "\n")

matrix.to_csv(out, sep="\t", mode="a")
