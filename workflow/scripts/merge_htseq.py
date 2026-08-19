#!/usr/bin/env python
"""Merge per-sample HTSeq-count outputs into a gene x sample matrix (drops __ rows).

kind = dedup  -> all samples with has_umi=true, from either assay branch
kind = 3prime -> legacy full-length 3'-restricted output (temporary compatibility)
Usage: merge_htseq.py <samples.tsv> <results_dir> <kind> <out.tsv>
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

    if kind == "dedup":
        if str(r["has_umi"]).strip().lower() not in _TRUE:
            continue
        if assay == "full_length_pe":
            branch = "full_length"
        elif assay == "quantseq_3prime_se":
            branch = "quantseq"
        else:
            continue
        suffix = ".dedup_htseq_counts.tsv"
    elif kind == "3prime":
        if assay != "full_length_pe":
            continue
        branch, suffix = "full_length", ".3prime_htseq_counts.tsv"
    else:
        sys.exit("unknown kind: %s" % kind)

    f = os.path.join(results, branch, s, s + suffix)
    d = pd.read_csv(f, sep="\t", header=None, names=["gene_id", "count"])
    d = d[~d["gene_id"].astype(str).str.startswith("__")]
    series[s] = d.set_index("gene_id")["count"].astype(int)

mat = pd.DataFrame(series).fillna(0).astype(int)
mat.index.name = "gene_id"
mat.to_csv(out, sep="\t")
