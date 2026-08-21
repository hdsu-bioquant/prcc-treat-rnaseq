#!/usr/bin/env python3
"""Concatenate stable per-library QC rows in sample-sheet order."""

import os
import sys
import pandas as pd

samples_f, results, out = sys.argv[1:4]
samples = pd.read_csv(samples_f, sep="\t", dtype=str, keep_default_na=False)
rows = []
for library in samples["library_id"]:
    path = os.path.join(results, "libraries", library, "qc_metrics.tsv")
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if len(table) != 1:
        raise SystemExit(f"Expected exactly one QC row in {path}")
    rows.append(table)
merged = pd.concat(rows, ignore_index=True)
os.makedirs(os.path.dirname(out), exist_ok=True)
merged.to_csv(out, sep="\t", index=False)
