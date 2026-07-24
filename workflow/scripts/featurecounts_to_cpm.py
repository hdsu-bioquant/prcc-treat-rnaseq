#!/usr/bin/env python
"""featureCounts table -> gene-level raw counts + CPM (QuantSeq 3' branch).

QuantSeq 3' tag counts are ~length-independent, so length-normalized metrics
(TPM/FPKM) are invalid; CPM is the correct within-sample normalization.
Ref (3' tag chemistry / no length normalization): Lexogen QuantSeq FWD documentation.
Usage: featurecounts_to_cpm.py <featureCounts.txt> <out.tsv>
Output columns: gene_id  length  count  cpm
"""
import sys
import pandas as pd

fc, out = sys.argv[1], sys.argv[2]
# featureCounts: first line is a '# Program...' comment, second line is the header.
df = pd.read_csv(fc, sep="\t", comment="#")
count_col = df.columns[-1]                       # last column = the BAM's counts
res = pd.DataFrame({
    "gene_id": df["Geneid"],
    "length":  df["Length"],
    "count":   df[count_col].astype(float),
})
total = res["count"].sum()
res["cpm"] = (res["count"] / total * 1e6) if total > 0 else 0.0
res.to_csv(out, sep="\t", index=False)
