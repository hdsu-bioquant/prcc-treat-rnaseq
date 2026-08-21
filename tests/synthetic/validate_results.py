#!/usr/bin/env python3
"""Validate the synthetic smoke test against the production output contract."""

from pathlib import Path
import hashlib
import sys
import pandas as pd

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
EXPRESSION_COLUMNS = [
    "gene_id", "gene_name", "gene_type", "gene_length",
    "unstranded", "stranded_first", "stranded_second",
    "fpkm", "fpkm_uq", "tpm", "umi_molecule_count",
]


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


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_checksum_file(path, base):
    if not path.is_file():
        fail(f"missing checksum file {path.relative_to(ROOT)}")
    checked = 0
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        digest, rel = line.split(None, 1)
        rel = rel.strip()
        target = base / rel
        if not target.is_file():
            fail(f"checksum target is missing: {target.relative_to(ROOT)}")
        if sha256(target) != digest:
            fail(f"checksum mismatch: {target.relative_to(ROOT)}")
        checked += 1
    if checked == 0:
        fail(f"checksum file is empty: {path.relative_to(ROOT)}")
    return checked


# Production-style top-level contract.
for directory in (RESULTS, RESTRICTED, INTERMEDIATE):
    if not directory.is_dir():
        fail(f"missing output directory {directory.relative_to(ROOT)}")
pass_("results/restricted/intermediate output structure")

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
    elif not pd.isna(row["umi_molecules_assigned"]):
        fail(f"{library} non-UMI QC should have NA umi_molecules_assigned")
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
# Validate stable identifiers generated by the custom-content table. Display text
# can be escaped or reformatted by different MultiQC versions.
if "prcc_rnaseq_qc" not in html_text or "prcc_rnaseq_qc_table" not in html_text:
    fail("library-level pRCC-RNA-Seq QC summary is missing from MultiQC report")
for label in ("Sample ID", "STAR input", "Unique mapped %", "UMI molecules"):
    if label not in html_text:
        fail(f"human-facing Library QC Summary label is missing from MultiQC report: {label}")
for software_name in software["software"].tolist():
    if str(software_name) not in html_text:
        fail(f"software version is missing from MultiQC report: {software_name}")
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

# Portable sample metadata must not contain raw FASTQ paths.
libraries = pd.read_csv(RESULTS / "run" / "libraries.tsv", sep="\t")
if "fq1" in libraries.columns or "fq2" in libraries.columns:
    fail("portable results/run/libraries.tsv unexpectedly contains FASTQ path columns")
pass_("portable library metadata omits raw FASTQ paths")

# Package integrity and deterministic validation manifests are self-consistent.
package_n = verify_checksum_file(RESULTS / "run" / "checksums.sha256", RESULTS)
validation_n = verify_checksum_file(RESULTS / "run" / "validation_checksums.sha256", RESULTS)
pass_(f"package checksums ({package_n} files)")
pass_(f"validation checksums ({validation_n} deterministic files)")

print("\nSynthetic smoke test PASSED.")
