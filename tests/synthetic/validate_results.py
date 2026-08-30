#!/usr/bin/env python3
"""Validate the synthetic smoke test against the production output contract."""

from pathlib import Path
from html import unescape
import argparse
import gzip
import hashlib
import re
import sys
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[2]
TEST = ROOT / "tests" / "synthetic"
OUTPUT = TEST / "output"
RESULTS = OUTPUT / "results"
RESTRICTED = OUTPUT / "restricted"
INTERMEDIATE = OUTPUT / "intermediate"
EXPECTED = TEST / "expected"

LIBRARIES = ["FL_noUMI", "FL_UMI", "QS_noUMI", "QS_UMI"]
UMI_LIBRARIES = ["FL_UMI", "QS_UMI"]
FL_LIBRARIES = ["FL_noUMI", "FL_UMI"]
QS_LIBRARIES = ["QS_noUMI", "QS_UMI"]
UMI_FIXTURE_SPECS = {
    "FL_UMI": {"pattern": "NNNNNNNN", "location": "read2_start", "discard_bases": 0},
    "QS_UMI": {"pattern": "NNNNNN", "location": "read1_start", "discard_bases": 4},
}
EXPRESSION_COLUMNS = [
    "gene_id", "gene_name", "gene_type", "gene_length",
    "unstranded", "stranded_first", "stranded_second",
    "fpkm", "fpkm_uq", "tpm", "umi_molecule_count",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Validate the pRCC-TREAT deterministic synthetic workflow outputs."
    )
    parser.add_argument(
        "--skip-frozen-baseline",
        action="store_true",
        help=(
            "validate fixture semantics, biological expectations, output contract and generated "
            "checksums, but do not compare results/run/validation_checksums.sha256 with the "
            "maintainer-frozen baseline. Intended only while deliberately changing the fixture/"
            "workflow before a reviewed re-freeze."
        ),
    )
    return parser.parse_args()


ARGS = parse_args()


def fail(msg):
    print(f"FAIL: {msg}", file=sys.stderr)
    raise SystemExit(1)


def pass_(msg):
    print(f"PASS: {msg}")


def read_expected(path, columns):
    df = pd.read_csv(path, sep="\t").set_index("gene_id")
    return df[columns].astype(int).sort_index()


def assert_equal(label, observed, expected):
    observed = observed.sort_index().sort_index(axis=1)
    expected = expected.sort_index().sort_index(axis=1)
    try:
        pd.testing.assert_frame_equal(observed, expected, check_dtype=False)
    except AssertionError as exc:
        print(f"FAIL: {label} differs from expected", file=sys.stderr)
        print(exc, file=sys.stderr)
        print("\nObserved:\n", observed, file=sys.stderr)
        print("\nExpected:\n", expected, file=sys.stderr)
        raise SystemExit(1)
    pass_(label)



def read_fastq_records(path):
    """Return FASTQ records as (header_without_at, sequence, quality)."""
    if not path.is_file():
        fail(f"missing FASTQ {path.relative_to(ROOT)}")

    records = []
    with gzip.open(path, "rt") as fh:
        line_no = 0
        while True:
            header = fh.readline()
            if not header:
                break
            line_no += 1
            sequence = fh.readline()
            plus = fh.readline()
            quality = fh.readline()
            if not sequence or not plus or not quality:
                fail(f"truncated FASTQ record beginning at line {line_no} in {path.relative_to(ROOT)}")
            line_no += 3

            header = header.rstrip("\n")
            sequence = sequence.rstrip("\n")
            plus = plus.rstrip("\n")
            quality = quality.rstrip("\n")
            if not header.startswith("@") or not plus.startswith("+"):
                fail(f"malformed FASTQ record beginning at line {line_no - 3} in {path.relative_to(ROOT)}")
            if len(sequence) != len(quality):
                fail(f"sequence/quality length mismatch in {path.relative_to(ROOT)} at {header}")
            records.append((header[1:], sequence, quality))
    return records


def raw_read_id(header):
    token = header.split()[0]
    return re.sub(r"/[12]$", "", token)


def validate_umi_fixture_definition():
    """Protect the intended synthetic UMI architectures before workflow execution."""
    samples = pd.read_csv(TEST / "samples.tsv", sep="\t", dtype=str, keep_default_na=False)
    required = {"library_id", "has_umi", "umi_pattern", "umi_location", "umi_discard_bases"}
    missing = sorted(required - set(samples.columns))
    if missing:
        fail("synthetic samples.tsv is missing UMI fixture column(s): " + ", ".join(missing))
    samples = samples.set_index("library_id")

    if set(samples.index) != set(LIBRARIES):
        fail(f"synthetic samples.tsv libraries differ from expected: {list(samples.index)}")

    for library, spec in UMI_FIXTURE_SPECS.items():
        row = samples.loc[library]
        if row["has_umi"].lower() != "true":
            fail(f"{library} must have has_umi=true in the synthetic fixture")
        observed = {
            "pattern": row["umi_pattern"],
            "location": row["umi_location"],
            "discard_bases": int(row["umi_discard_bases"]),
        }
        if observed != spec:
            fail(f"{library} UMI fixture spec differs (observed={observed}, expected={spec})")

    for library in sorted(set(LIBRARIES) - set(UMI_LIBRARIES)):
        row = samples.loc[library]
        if row["has_umi"].lower() != "false":
            fail(f"{library} must have has_umi=false in the synthetic fixture")
        if any(row[col] not in {"", "-"} for col in ("umi_pattern", "umi_location", "umi_discard_bases")):
            fail(f"{library} non-UMI fixture row contains substantive UMI metadata")

    manifest = pd.read_csv(EXPECTED / "read_manifest.tsv", sep="\t", dtype=str).set_index("read_id")

    fl_r1 = read_fastq_records(TEST / "data" / "FL_UMI_R1.fastq.gz")
    fl_r2 = read_fastq_records(TEST / "data" / "FL_UMI_R2.fastq.gz")
    if len(fl_r1) != len(fl_r2):
        fail("FL_UMI raw R1/R2 record counts differ")
    for r1, r2 in zip(fl_r1, fl_r2):
        rid1 = raw_read_id(r1[0])
        rid2 = raw_read_id(r2[0])
        if rid1 != rid2:
            fail(f"FL_UMI paired raw read IDs differ: {r1[0]} vs {r2[0]}")
        expected_umi = manifest.loc[rid1, "umi"]
        if r2[1][:8] != expected_umi:
            fail(f"FL_UMI R2 does not start with expected 8-nt UMI for {rid1}")
        if len(r1[1]) != 94 or len(r2[1]) != 108:
            fail(f"FL_UMI raw read lengths do not encode the intended R2-start UMI for {rid1}")

    qs = read_fastq_records(TEST / "data" / "QS_UMI_R1.fastq.gz")
    for header, sequence, _quality in qs:
        rid = raw_read_id(header)
        expected_umi = manifest.loc[rid, "umi"]
        if sequence[:6] != expected_umi:
            fail(f"QS_UMI R1 does not start with expected 6-nt UMI for {rid}")
        if sequence[6:10] != "TATA":
            fail(f"QS_UMI R1 is missing the 4-nt TATA spacer at bases 7-10 for {rid}")
        if len(sequence) != 79:
            fail(f"QS_UMI raw read length is not 79 nt for {rid}")

    pass_("synthetic UMI fixture definitions (R2/8+0 and R1/6+4 TATA)")


def assert_umi_in_processed_header(header, umi, library):
    # Raw synthetic names deliberately do not contain the UMI sequence. UMI-tools
    # must therefore add it to the processed name for downstream deduplication.
    token = header.split()[0]
    if umi not in token:
        fail(f"{library} processed read name does not contain extracted UMI {umi}: {header}")


def validate_umi_extraction_intermediates():
    """Verify exact raw->UMI-extracted transformations, independent of STAR soft clipping."""
    manifest = pd.read_csv(EXPECTED / "read_manifest.tsv", sep="\t", dtype=str).set_index("read_id")

    raw_fl_r1 = read_fastq_records(TEST / "data" / "FL_UMI_R1.fastq.gz")
    raw_fl_r2 = read_fastq_records(TEST / "data" / "FL_UMI_R2.fastq.gz")
    out_fl_r1 = read_fastq_records(INTERMEDIATE / "libraries" / "FL_UMI" / "preprocess" / "R1.umi.fastq.gz")
    out_fl_r2 = read_fastq_records(INTERMEDIATE / "libraries" / "FL_UMI" / "preprocess" / "R2.umi.fastq.gz")
    if not (len(raw_fl_r1) == len(raw_fl_r2) == len(out_fl_r1) == len(out_fl_r2)):
        fail("FL_UMI record count changed during UMI extraction")
    for raw1, raw2, out1, out2 in zip(raw_fl_r1, raw_fl_r2, out_fl_r1, out_fl_r2):
        rid = raw_read_id(raw1[0])
        if raw_read_id(raw2[0]) != rid:
            fail(f"FL_UMI raw pair IDs differ for {rid}")
        umi = manifest.loc[rid, "umi"]
        if out1[1] != raw1[1] or out1[2] != raw1[2]:
            fail(f"FL_UMI R1 was unexpectedly altered during R2-start UMI extraction for {rid}")
        if out2[1] != raw2[1][8:] or out2[2] != raw2[2][8:]:
            fail(f"FL_UMI R2 did not remove exactly the 8-nt UMI for {rid}")
        assert_umi_in_processed_header(out1[0], umi, "FL_UMI")
        assert_umi_in_processed_header(out2[0], umi, "FL_UMI")

    raw_qs = read_fastq_records(TEST / "data" / "QS_UMI_R1.fastq.gz")
    out_qs = read_fastq_records(INTERMEDIATE / "libraries" / "QS_UMI" / "preprocess" / "R1.umi.fastq.gz")
    if len(raw_qs) != len(out_qs):
        fail("QS_UMI record count changed during UMI extraction")
    for raw, out in zip(raw_qs, out_qs):
        rid = raw_read_id(raw[0])
        umi = manifest.loc[rid, "umi"]
        if raw[1][6:10] != "TATA":
            fail(f"QS_UMI raw fixture spacer unexpectedly differs from TATA for {rid}")
        if out[1] != raw[1][10:] or out[2] != raw[2][10:]:
            fail(f"QS_UMI did not remove exactly 6 UMI + 4 discard bases for {rid}")
        assert_umi_in_processed_header(out[0], umi, "QS_UMI")

    pass_("exact UMI extraction intermediates (including QuantSeq spacer removal)")


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_checksum_file(path):
    """Read a sha256sum-style manifest as {relative_path: digest}."""
    if not path.is_file():
        fail(f"missing checksum file {path.relative_to(ROOT)}")

    entries = {}
    for line_no, line in enumerate(path.read_text().splitlines(), start=1):
        if not line.strip():
            continue
        try:
            digest, rel = line.split(None, 1)
        except ValueError:
            fail(f"malformed checksum line {line_no} in {path.relative_to(ROOT)}")
        rel = rel.strip()
        if not re.fullmatch(r"[0-9a-fA-F]{64}", digest):
            fail(f"invalid SHA256 digest on line {line_no} in {path.relative_to(ROOT)}")
        if rel in entries:
            fail(f"duplicate checksum target {rel} in {path.relative_to(ROOT)}")
        entries[rel] = digest.lower()

    if not entries:
        fail(f"checksum file is empty: {path.relative_to(ROOT)}")
    return entries


def verify_checksum_file(path, base):
    entries = parse_checksum_file(path)
    for rel, digest in entries.items():
        target = base / rel
        if not target.is_file():
            fail(f"checksum target is missing: {target.relative_to(ROOT)}")
        observed = sha256(target)
        if observed != digest:
            fail(
                f"checksum mismatch: {target.relative_to(ROOT)} "
                f"(observed={observed}, expected={digest})"
            )
    return len(entries)


def compare_checksum_files(observed_path, expected_path):
    """Compare generated validation hashes with the frozen reference baseline."""
    observed = parse_checksum_file(observed_path)
    expected = parse_checksum_file(expected_path)

    observed_paths = set(observed)
    expected_paths = set(expected)
    missing = sorted(expected_paths - observed_paths)
    unexpected = sorted(observed_paths - expected_paths)
    changed = sorted(
        rel for rel in observed_paths & expected_paths
        if observed[rel] != expected[rel]
    )

    if missing or unexpected or changed:
        print("FAIL: frozen validation checksum baseline differs", file=sys.stderr)
        if missing:
            print("  Missing deterministic files:", file=sys.stderr)
            for rel in missing:
                print(f"    {rel}", file=sys.stderr)
        if unexpected:
            print("  Unexpected deterministic files:", file=sys.stderr)
            for rel in unexpected:
                print(f"    {rel}", file=sys.stderr)
        if changed:
            print("  Changed deterministic files:", file=sys.stderr)
            for rel in changed:
                print(f"    {rel}", file=sys.stderr)
                print(f"      observed: {observed[rel]}", file=sys.stderr)
                print(f"      expected: {expected[rel]}", file=sys.stderr)
        raise SystemExit(1)

    return len(expected)


validate_umi_fixture_definition()

# Production-style top-level contract.
for directory in (RESULTS, RESTRICTED, INTERMEDIATE):
    if not directory.is_dir():
        fail(f"missing output directory {directory.relative_to(ROOT)}")
pass_("results/restricted/intermediate output structure")

# UMI extraction intermediates are retained by run_test.sh via --notemp so the
# regression test can detect incorrect 5-prime preprocessing before alignment.
validate_umi_extraction_intermediates()

# Exact raw-count matrix.
raw_path = RESULTS / "matrices" / "raw_gene_counts.tsv"
if not raw_path.exists():
    fail(f"missing {raw_path.relative_to(ROOT)}")
raw = pd.read_csv(raw_path, sep="\t").set_index("gene_id")
if list(raw.columns) != LIBRARIES:
    fail(f"raw matrix columns are {list(raw.columns)}, expected {LIBRARIES}")
raw = raw.astype(int)
exp_raw = read_expected(EXPECTED / "raw_gene_counts.tsv", LIBRARIES)
assert_equal("raw STAR gene-count matrix", raw, exp_raw)

# Exact UMI molecule matrix.
umi_path = RESULTS / "matrices" / "umi_molecule_counts.tsv"
if not umi_path.exists():
    fail(f"missing {umi_path.relative_to(ROOT)}")
umi = pd.read_csv(umi_path, sep="\t").set_index("gene_id")
if list(umi.columns) != UMI_LIBRARIES:
    fail(f"UMI matrix columns are {list(umi.columns)}, expected {UMI_LIBRARIES}")
umi = umi.astype(int)
exp_umi = read_expected(EXPECTED / "umi_dedup_gene_counts.tsv", UMI_LIBRARIES)
assert_equal("UMI molecule-count matrix", umi, exp_umi)

# Canonical per-library expression files use one schema for both assays.
expected_lengths = pd.read_csv(EXPECTED / "gene_lengths.tsv", sep="\t").set_index("gene_id")
for library in LIBRARIES:
    path = RESULTS / "libraries" / library / "gene_expression.tsv"
    if not path.is_file():
        fail(f"missing {path.relative_to(ROOT)}")
    expr = pd.read_csv(path, sep="\t", na_values=["NA"], keep_default_na=True)
    if list(expr.columns) != EXPRESSION_COLUMNS:
        fail(f"{library} gene_expression.tsv has unexpected columns: {list(expr.columns)}")
    expr = expr.set_index("gene_id")

    obs = expr[["unstranded"]].rename(columns={"unstranded": library}).astype(int)
    exp = exp_raw[[library]]
    assert_equal(f"{library} per-library unstranded counts", obs, exp)

    for gid in expected_lengths.index:
        if gid not in expr.index:
            fail(f"{library} expression table is missing {gid}")
        if int(expr.loc[gid, "gene_length"]) != int(expected_lengths.loc[gid, "gene_length"]):
            fail(f"{library} has wrong gene length for {gid}")
        if expr.loc[gid, "gene_name"] != expected_lengths.loc[gid, "gene_name"]:
            fail(f"{library} has wrong gene_name for {gid}")
        if expr.loc[gid, "gene_type"] != expected_lengths.loc[gid, "gene_type"]:
            fail(f"{library} has wrong gene_type for {gid}")

    if library in FL_LIBRARIES:
        if expr[["fpkm", "fpkm_uq", "tpm"]].isna().any().any():
            fail(f"{library} full-length normalized expression contains NA")
        if abs(float(expr["tpm"].sum()) - 1_000_000.0) > 1.0:
            fail(f"{library} TPM values do not sum to approximately 1e6")
    else:
        if not expr[["fpkm", "fpkm_uq", "tpm"]].isna().all().all():
            fail(f"{library} QuantSeq length-normalized expression should be NA")

    if library in UMI_LIBRARIES:
        observed = expr[["umi_molecule_count"]].rename(columns={"umi_molecule_count": library}).astype(int)
        assert_equal(f"{library} per-library UMI molecule counts", observed, exp_umi[[library]])
    elif not expr["umi_molecule_count"].isna().all():
        fail(f"{library} non-UMI expression table should have NA umi_molecule_count")
pass_("canonical gene_expression.tsv schema")

# Stable QC contract and selected exact fixture expectations.
summary = pd.read_csv(EXPECTED / "expected_summary.tsv", sep="\t").set_index("library_id")
run_qc_path = RESULTS / "qc" / "qc_metrics.tsv"
if not run_qc_path.is_file():
    fail(f"missing {run_qc_path.relative_to(ROOT)}")
run_qc = pd.read_csv(run_qc_path, sep="\t", na_values=["NA"], keep_default_na=True).set_index("library_id")
if list(run_qc.index) != LIBRARIES:
    fail(f"run QC libraries are {list(run_qc.index)}, expected {LIBRARIES}")
required_qc_columns = {
    "umi_length",
    "umi_location",
    "umi_discard_bases",
    "umi_extract_qc_records",
    "umi_extract_qc_retained_percent",
    "umi_extract_transform_match_percent",
    "umi_extract_tag_match_percent",
}
missing_qc_columns = sorted(required_qc_columns - set(run_qc.columns))
if missing_qc_columns:
    fail("run QC table is missing UMI extraction QC column(s): " + ", ".join(missing_qc_columns))
for library in LIBRARIES:
    per_lib = RESULTS / "libraries" / library / "qc_metrics.tsv"
    if not per_lib.is_file():
        fail(f"missing {per_lib.relative_to(ROOT)}")
    row = run_qc.loc[library]
    observed_input = int(row["star_input_records"])
    expected_input = int(summary.loc[library, "expected_after_trim"])
    if observed_input != expected_input:
        fail(
            f"{library} star_input_records differs from expected "
            f"(observed={observed_input}, expected={expected_input})"
        )
    if int(row["gene_assigned_unstranded"]) != int(summary.loc[library, "expected_raw_gene_assigned"]):
        fail(f"{library} gene_assigned_unstranded differs from expected")
    if library in UMI_LIBRARIES:
        if int(row["umi_molecules_assigned"]) != int(summary.loc[library, "expected_umi_molecules"]):
            fail(f"{library} umi_molecules_assigned differs from expected")
        spec = UMI_FIXTURE_SPECS[library]
        if int(row["umi_length"]) != len(spec["pattern"]):
            fail(f"{library} QC umi_length differs from fixture specification")
        if str(row["umi_location"]) != spec["location"]:
            fail(f"{library} QC umi_location differs from fixture specification")
        if int(row["umi_discard_bases"]) != spec["discard_bases"]:
            fail(f"{library} QC umi_discard_bases differs from fixture specification")
        expected_checked = 16 if library == "FL_UMI" else 17
        if int(row["umi_extract_qc_records"]) != expected_checked:
            fail(f"{library} UMI extraction QC checked-record count differs from fixture")
        for col in (
            "umi_extract_qc_retained_percent",
            "umi_extract_transform_match_percent",
            "umi_extract_tag_match_percent",
        ):
            if float(row[col]) != 100.0:
                fail(f"{library} {col} should be 100% in the deterministic fixture")
    else:
        if not pd.isna(row["umi_molecules_assigned"]):
            fail(f"{library} non-UMI QC should have NA umi_molecules_assigned")
        for col in required_qc_columns:
            if not pd.isna(row[col]):
                fail(f"{library} non-UMI QC should have NA {col}")
pass_("canonical per-library and run-level QC metrics")

# Software provenance is manifest-driven (one production rule, no per-tool probes)
# plus the actual Snakemake controller version from the run environment.
software_path = RESULTS / "run" / "software_versions.tsv"
if not software_path.is_file():
    fail("missing results/run/software_versions.tsv")
software = pd.read_csv(software_path, sep="\t").set_index("tool_id")
expected_tools = [
    "bbmap", "fastqc", "htseq", "multiqc", "py", "samtools", "star", "umitools",
    "snakemake",
]
if sorted(software.index.tolist()) != sorted(expected_tools):
    fail(
        "software_versions.tsv tool set differs from expected "
        f"(observed={sorted(software.index.tolist())}, expected={sorted(expected_tools)})"
    )
required_software_columns = {
    "software", "version", "version_source", "container_source", "resolved_container"
}
if not required_software_columns.issubset(software.columns):
    fail("software_versions.tsv is missing required provenance columns")
for tool in expected_tools:
    if not str(software.loc[tool, "version"]).strip():
        fail(f"software_versions.tsv has an empty version for {tool}")
if str(software.loc["umitools", "version"]) != "1.1.6":
    fail("synthetic qualification requires maintained UMI-tools version 1.1.6")
if str(software.loc["umitools", "container_source"]) != (
    "docker://quay.io/biocontainers/umi_tools:1.1.6--py39hbcbf7aa_0"
):
    fail("synthetic qualification UMI-tools container source differs from maintained 1.1.6 image")
for tool in expected_tools:
    if tool == "snakemake":
        if software.loc[tool, "version_source"] != "runtime":
            fail("Snakemake version should be recorded from the runtime controller environment")
        if str(software.loc[tool, "version"]).strip().lower() == "unknown":
            fail("Snakemake runtime version could not be determined")
    elif software.loc[tool, "version_source"] != "pipeline_manifest":
        fail(f"{tool} version should come from the pipeline software manifest")
pass_("software-version provenance")

# Human-facing MultiQC report.
multiqc = RESULTS / "qc" / "multiqc_report.html"
if not multiqc.is_file() or multiqc.stat().st_size == 0:
    fail("MultiQC report is missing or empty")
if not (RESTRICTED / "qc" / "multiqc_data").is_dir():
    fail("restricted MultiQC data directory is missing")
html_text = multiqc.read_text(errors="replace")
# Validate stable identifiers generated by the custom-content table. MultiQC can
# HTML-escape or wrap visible labels, so normalise rendered text before checking
# presentation strings instead of relying on raw HTML byte-for-byte text.
if "prcc_treat_rnaseq_qc" not in html_text or "prcc_treat_rnaseq_qc_table" not in html_text:
    fail("library-level pRCC-TREAT-RNA-Seq QC summary is missing from MultiQC report")
rendered_text = unescape(re.sub(r"<[^>]+>", " ", html_text))
rendered_text = " ".join(rendered_text.split())
for label in (
    "Library ID", "Sample ID", "STAR input", "Unique mapped %",
    "UMI design", "UMI transform %", "UMI molecules",
):
    if label not in rendered_text:
        fail(f"human-facing Library QC Summary label is missing from MultiQC report: {label}")

# MultiQC also exports its software-version data under multiqc_data. The exact
# TSV layout is a MultiQC implementation detail and is not part of our pipeline
# contract, so do not parse it as a fixed two-column table. MultiQC 1.21 can also
# normalise version strings for display (for example 39.06 -> 39.6 and
# 2.7.5c -> 2.7.5rc0). Exact provenance is already validated above from the
# pipeline-controlled results/run/software_versions.tsv. Here we only verify that
# the MultiQC export exists and that every declared software name reached both
# the machine-readable MultiQC export and the rendered report.
version_export = RESTRICTED / "qc" / "multiqc_data" / "multiqc_software_versions.txt"
if not version_export.is_file() or version_export.stat().st_size == 0:
    fail("MultiQC software-version export is missing")
version_export_text = version_export.read_text(errors="replace")

for _, row in software.iterrows():
    software_name = str(row["software"])
    if software_name not in version_export_text:
        fail(f"MultiQC software-version export is missing software entry: {software_name}")
    if software_name not in rendered_text:
        fail(f"software name is missing from MultiQC report: {software_name}")
pass_("MultiQC report + library QC summary + complete software versions + site-retained MultiQC data")

# Restricted site-retained alignment/QC products.
for library in LIBRARIES:
    bam = RESTRICTED / "libraries" / library / "alignments" / "genomic.sorted.bam"
    bai = RESTRICTED / "libraries" / library / "alignments" / "genomic.sorted.bam.bai"
    fastqc = RESTRICTED / "libraries" / library / "qc" / "fastqc"
    if not bam.is_file() or bam.stat().st_size == 0 or not bai.is_file():
        fail(f"restricted genomic BAM/index missing for {library}")
    if not fastqc.is_dir() or not list(fastqc.glob("*_fastqc.html")):
        fail(f"restricted FastQC output missing for {library}")
for library in UMI_LIBRARIES:
    bam = RESTRICTED / "libraries" / library / "alignments" / "umi_dedup.bam"
    if not bam.is_file() or bam.stat().st_size == 0:
        fail(f"restricted UMI-deduplicated BAM missing for {library}")
pass_("restricted alignment and detailed QC products")

# Portable and restricted run metadata.
portable_run = [
    "libraries.tsv", "config.yaml", "provenance.yaml", "software_versions.tsv",
    "references.tsv", "manifest.tsv", "checksums.sha256", "validation_checksums.sha256",
]
for name in portable_run:
    if not (RESULTS / "run" / name).is_file():
        fail(f"missing results/run/{name}")
for name in ("libraries.original.tsv", "config.effective.yaml"):
    if not (RESTRICTED / "run" / name).is_file():
        fail(f"missing restricted/run/{name}")
pass_("portable and restricted run metadata")

# Portable pipeline identity must come from maintainer-owned release metadata.
with (RESULTS / "run" / "config.yaml").open() as fh:
    portable_config = yaml.safe_load(fh)
with (ROOT / "workflow" / "release.yaml").open() as fh:
    release_metadata = yaml.safe_load(fh)
expected_pipeline = {
    "name": release_metadata["pipeline_name"],
    "release": release_metadata["pipeline_release"],
    "output_contract": int(release_metadata["output_contract"]),
}
if portable_config.get("pipeline") != expected_pipeline:
    fail(
        "portable config pipeline identity differs from maintainer-owned workflow/release.yaml: "
        f"{portable_config.get('pipeline')}"
    )
pass_("maintainer-owned pipeline release metadata")

# Portable sample metadata must not contain raw FASTQ paths.
libraries = pd.read_csv(RESULTS / "run" / "libraries.tsv", sep="\t")
if "fq1" in libraries.columns or "fq2" in libraries.columns:
    fail("portable results/run/libraries.tsv unexpectedly contains FASTQ path columns")
pass_("portable library metadata omits raw FASTQ paths")

# Package integrity and deterministic validation manifests are self-consistent.
package_n = verify_checksum_file(RESULTS / "run" / "checksums.sha256", RESULTS)
validation_path = RESULTS / "run" / "validation_checksums.sha256"
validation_n = verify_checksum_file(validation_path, RESULTS)
pass_(f"package checksums ({package_n} files)")
pass_(f"validation checksums ({validation_n} deterministic files)")

# Cross-installation qualification: by default the deterministic output subset must
# match the maintainer-approved frozen reference. During an intentional refactor,
# --skip-frozen-baseline permits all other validation to complete before re-freezing.
if ARGS.skip_frozen_baseline:
    pass_("frozen cross-installation validation baseline comparison intentionally skipped")
else:
    frozen_validation = EXPECTED / "validation_checksums.sha256"
    frozen_n = compare_checksum_files(validation_path, frozen_validation)
    if frozen_n != validation_n:
        fail(
            "frozen validation checksum entry count differs from generated validation set "
            f"(generated={validation_n}, frozen={frozen_n})"
        )
    pass_(f"frozen cross-installation validation baseline ({frozen_n} files)")

print("\nSynthetic smoke test PASSED.")
