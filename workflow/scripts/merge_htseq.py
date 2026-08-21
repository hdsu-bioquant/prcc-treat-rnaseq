#!/usr/bin/env python3
"""Merge per-library UMI-deduplicated HTSeq counts into a gene x library matrix.

Drops HTSeq summary rows beginning with ``__`` and includes every library with
``has_umi=true`` regardless of assay.

Usage: merge_htseq.py <samples.tsv> <results_dir> <kind> <out.tsv>
Currently supported kind: dedup
"""

import os
import sys

import pandas as pd


samples_f, results, kind, out = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
if kind != "dedup":
    sys.exit("unknown kind: %s (supported: dedup)" % kind)

samples = pd.read_csv(samples_f, sep="\t", dtype=str, keep_default_na=False).set_index("library_id", drop=False)

series = {}
for library_id, row in samples.iterrows():
    if row["has_umi"] != "true":
        continue

    assay = row["assay"]
    branch = "full_length" if assay == "full_length" else "quantseq"
    path = os.path.join(results, branch, library_id, library_id + ".dedup_htseq_counts.tsv")
    counts = pd.read_csv(path, sep="\t", header=None, names=["gene_id", "count"])
    counts = counts[~counts["gene_id"].astype(str).str.startswith("__")]
    series[library_id] = counts.set_index("gene_id")["count"].astype(int)

matrix = pd.DataFrame(series).fillna(0).astype(int)
matrix.index.name = "gene_id"
matrix.to_csv(out, sep="\t")
