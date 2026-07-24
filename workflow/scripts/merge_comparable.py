#!/usr/bin/env python
"""Cross-assay comparable CPM matrix (gene x sample).

Pulls the CPM column from each sample's length-bias-matched output:
  full-length -> {s}.3prime_gene_counts_cpm.tsv  (3'-window CPM, emulating 3'-tag)
  QuantSeq    -> {s}.gene_counts_cpm.tsv          (already 3'/length-independent CPM)
The `assay` annotation row is carried so downstream comparisons stay assay-aware.
Usage: merge_comparable.py <samples.tsv> <results_dir> <out.tsv>
"""
import sys, os
import pandas as pd

samples_f, results, out = sys.argv[1], sys.argv[2], sys.argv[3]
samples = pd.read_csv(samples_f, sep="\t", dtype=str).fillna("-")

series, assays = {}, {}
for _, r in samples.iterrows():
    s, assay = r["sample"], r["assay"]
    if assay == "full_length_pe":
        f = os.path.join(results, "full_length", s, s + ".3prime_gene_counts_cpm.tsv")
    elif assay == "quantseq_3prime_se":
        f = os.path.join(results, "quantseq", s, s + ".gene_counts_cpm.tsv")
    else:
        continue
    d = pd.read_csv(f, sep="\t")
    series[s] = d.set_index("gene_id")["cpm"]
    assays[s] = assay

mat = pd.DataFrame(series).fillna(0)
mat.index.name = "gene_id"
with open(out, "w") as o:
    o.write("# assay\t" + "\t".join(assays[c] for c in mat.columns) + "\n")
    o.write("# values = CPM on a common 3'-comparable gene-count basis; compare WITHIN gene across samples,\n")
    o.write("# and prefer relative/fold-change comparisons across assays (see docs/Cross_Assay_Comparability).\n")
mat.to_csv(out, sep="\t", mode="a")
