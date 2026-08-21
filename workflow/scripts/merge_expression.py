#!/usr/bin/env python3
"""Merge canonical per-library gene-expression tables into a gene x library matrix.

Usage:
  merge_expression.py <samples.tsv> <results_dir> <metric> <out.tsv>

Supported metrics:
  unstranded         all libraries
  umi_molecule_count only has_umi=true libraries
"""

import os
import sys
import pandas as pd

samples_f, results, metric, out = sys.argv[1:5]
if metric not in {"unstranded", "umi_molecule_count"}:
    raise SystemExit(f"Unsupported metric: {metric}")

samples = pd.read_csv(samples_f, sep="\t", dtype=str, keep_default_na=False)
if metric == "umi_molecule_count":
    samples = samples[samples["has_umi"].str.lower() == "true"]

series = {}
for _, row in samples.iterrows():
    library = row["library_id"]
    path = os.path.join(results, "libraries", library, "gene_expression.tsv")
    table = pd.read_csv(path, sep="\t", na_values=["NA"], keep_default_na=True)
    if metric not in table.columns:
        raise SystemExit(f"{path} is missing metric column '{metric}'")
    values = table.set_index("gene_id")[metric]
    if metric == "umi_molecule_count" and values.isna().any():
        raise SystemExit(f"{path} contains NA UMI molecule counts for a UMI-bearing library")
    series[library] = values.astype(int)

matrix = pd.DataFrame(series).fillna(0).astype(int)
matrix.index.name = "gene_id"
os.makedirs(os.path.dirname(out), exist_ok=True)
matrix.to_csv(out, sep="\t")
