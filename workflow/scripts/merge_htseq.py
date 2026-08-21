#!/usr/bin/env python
"""Merge per-sample UMI-deduplicated HTSeq-count outputs into a gene x sample matrix.

Drops HTSeq summary rows beginning with ``__`` and includes every library with
``has_umi=true`` regardless of assay.
Usage: merge_htseq.py <samples.tsv> <results_dir> <kind> <out.tsv>
Currently supported kind: dedup
"""
import sys, os
import pandas as pd

samples_f, results, kind, out = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
samples = pd.read_csv(samples_f, sep="\t", dtype=str).fillna("-")
_TRUE = {"true", "1", "yes", "y"}

# Backwards compatibility for old sheets during the transition. The old pipeline
# had one global QuantSeq UMI toggle, so if has_umi is absent and this rule is
# requested, QuantSeq rows are the legacy UMI-bearing rows.
if "has_umi" not in samples.columns:
    samples["has_umi"] = "false"
    samples.loc[samples["assay"] == "quantseq_3prime_se", "has_umi"] = "true"

series = {}
for _, r in samples.iterrows():
    s = r["sample"]
    assay = r["assay"]

    if kind != "dedup":
        sys.exit("unknown kind: %s (supported: dedup)" % kind)

    if str(r["has_umi"]).strip().lower() not in _TRUE:
        continue
    if assay == "full_length_pe":
        branch = "full_length"
    elif assay == "quantseq_3prime_se":
        branch = "quantseq"
    else:
        continue
    suffix = ".dedup_htseq_counts.tsv"

    f = os.path.join(results, branch, s, s + suffix)
    d = pd.read_csv(f, sep="\t", header=None, names=["gene_id", "count"])
    d = d[~d["gene_id"].astype(str).str.startswith("__")]
    series[s] = d.set_index("gene_id")["count"].astype(int)

mat = pd.DataFrame(series).fillna(0).astype(int)
mat.index.name = "gene_id"
mat.to_csv(out, sep="\t")
