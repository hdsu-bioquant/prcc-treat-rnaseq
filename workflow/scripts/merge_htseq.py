#!/usr/bin/env python
"""Merge per-sample HTSeq-count outputs into a gene x sample matrix (drops __ rows).

Used for the two SECONDARY matrices:
  kind = dedup   -> quantseq/{s}/{s}.dedup_htseq_counts.tsv     (QuantSeq UMI-deduplicated)
  kind = 3prime  -> full_length/{s}/{s}.3prime_htseq_counts.tsv (Branch-A 3'-restricted)
Usage: merge_htseq.py <samples.tsv> <results_dir> <kind> <out.tsv>
"""
import sys, os
import pandas as pd

samples_f, results, kind, out = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
samples = pd.read_csv(samples_f, sep="\t", dtype=str).fillna("-")

if kind == "dedup":
    sel, branch, suffix = "quantseq_3prime_se", "quantseq", ".dedup_htseq_counts.tsv"
elif kind == "3prime":
    sel, branch, suffix = "full_length_pe", "full_length", ".3prime_htseq_counts.tsv"
else:
    sys.exit("unknown kind: %s" % kind)

series = {}
for _, r in samples.iterrows():
    if r["assay"] != sel:
        continue
    s = r["sample"]
    f = os.path.join(results, branch, s, s + suffix)
    d = pd.read_csv(f, sep="\t", header=None, names=["gene_id", "count"])
    d = d[~d["gene_id"].astype(str).str.startswith("__")]   # drop __no_feature/__ambiguous/...
    series[s] = d.set_index("gene_id")["count"].astype(int)

mat = pd.DataFrame(series).fillna(0).astype(int)
mat.index.name = "gene_id"
mat.to_csv(out, sep="\t")
