#!/usr/bin/env python3
"""Validate the synthetic smoke-test's canonical numerical outputs.

Run after Snakemake finishes:
    python tests/synthetic/validate_results.py

The check intentionally compares parsed numerical tables, not whole-file hashes.
"""
from pathlib import Path
import sys
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
TEST = ROOT / "tests" / "synthetic"
RESULTS = TEST / "results"
EXPECTED = TEST / "expected"


def fail(msg):
    print(f"FAIL: {msg}", file=sys.stderr)
    raise SystemExit(1)


def read_expected(path, sample_cols):
    df = pd.read_csv(path, sep="\t").set_index("gene_id")
    return df[sample_cols].astype(int).sort_index()


def assert_equal(label, observed, expected):
    observed = observed.sort_index().sort_index(axis=1)
    expected = expected.sort_index().sort_index(axis=1)
    try:
        pd.testing.assert_frame_equal(observed, expected, check_dtype=False)
    except AssertionError as e:
        print(f"FAIL: {label} differs from expected", file=sys.stderr)
        print(e, file=sys.stderr)
        print("\nObserved:\n", observed, file=sys.stderr)
        print("\nExpected:\n", expected, file=sys.stderr)
        raise SystemExit(1)
    print(f"PASS: {label}")


raw_path = RESULTS / "matrix" / "gene_counts_matrix.tsv"
if not raw_path.exists():
    fail(f"missing {raw_path.relative_to(ROOT)}")
raw = pd.read_csv(raw_path, sep="\t", comment="#").set_index("gene_id")
raw_samples = ["FL_noUMI", "FL_UMI", "QS_noUMI", "QS_UMI"]
missing = [s for s in raw_samples if s not in raw.columns]
if missing:
    fail(f"raw matrix is missing samples: {', '.join(missing)}")
raw = raw[raw_samples].astype(int)
exp_raw = read_expected(EXPECTED / "raw_gene_counts.tsv", raw_samples)
assert_equal("raw STAR gene counts", raw, exp_raw)

umi_path = RESULTS / "matrix" / "umi_dedup_matrix.tsv"
if not umi_path.exists():
    fail(f"missing {umi_path.relative_to(ROOT)}")
umi = pd.read_csv(umi_path, sep="\t").set_index("gene_id")
umi_samples = ["FL_UMI", "QS_UMI"]
if set(umi.columns) != set(umi_samples):
    fail(f"UMI matrix columns are {list(umi.columns)}, expected exactly {umi_samples}")
umi = umi[umi_samples].astype(int)
exp_umi = read_expected(EXPECTED / "umi_dedup_gene_counts.tsv", umi_samples)
assert_equal("UMI-deduplicated molecule counts", umi, exp_umi)

multiqc = RESULTS / "qc" / "multiqc_report.html"
if not multiqc.exists() or multiqc.stat().st_size == 0:
    fail("MultiQC report is missing or empty")
print("PASS: MultiQC report exists")

print("\nSynthetic smoke test PASSED.")
