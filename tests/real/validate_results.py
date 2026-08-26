#!/usr/bin/env python3
"""Validate a completed realistic public-data qualification run.

The qualification is intentionally documentation-led: users manually copy and edit the
maintained tests/real/config.yaml + samples.tsv and configure/reuse a site execution profile.
This validator checks that the resulting run still matches the pinned qualification definition
and that the generated portable/restricted outputs are coherent.

During initial maintainer qualification, use --skip-frozen-baseline until two clean runs have
reproduced the same results/run/validation_checksums.sha256 and the baseline is deliberately
frozen under tests/real/expected/.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import re
import sys
from collections import Counter
from html import unescape
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[2]
TEST = ROOT / "tests" / "real"
EXPECTED = TEST / "expected"
LIBRARIES = ["REAL_FL", "REAL_QS_UMI"]
UMI_LIBRARIES = ["REAL_QS_UMI"]
EXPRESSION_COLUMNS = [
    "gene_id", "gene_name", "gene_type", "gene_length",
    "unstranded", "stranded_first", "stranded_second",
    "fpkm", "fpkm_uq", "tpm", "umi_molecule_count",
]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--run-dir",
        required=True,
        help="Qualification run directory containing the copied config.yaml and samples.tsv",
    )
    parser.add_argument(
        "--skip-frozen-baseline",
        action="store_true",
        help="Skip only comparison with tests/real/expected/validation_checksums.sha256",
    )
    return parser.parse_args()


ARGS = parse_args()
RUN_DIR = Path(ARGS.run_dir).expanduser().resolve()
CONFIG_PATH = RUN_DIR / "config.yaml"
SAMPLES_PATH = RUN_DIR / "samples.tsv"


def fail(message):
    print(f"FAIL: {message}", file=sys.stderr)
    raise SystemExit(1)


def pass_(message):
    print(f"PASS: {message}")


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def md5(path):
    h = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_checksum_file(path):
    if not path.is_file():
        fail(f"missing checksum file: {path}")
    entries = {}
    for lineno, line in enumerate(path.read_text().splitlines(), 1):
        if not line.strip():
            continue
        try:
            digest, rel = line.split(None, 1)
        except ValueError:
            fail(f"malformed checksum line {lineno} in {path}")
        rel = rel.strip()
        if not re.fullmatch(r"[0-9a-fA-F]{64}", digest):
            fail(f"invalid SHA256 on line {lineno} in {path}")
        if rel in entries:
            fail(f"duplicate checksum target {rel} in {path}")
        entries[rel] = digest.lower()
    if not entries:
        fail(f"checksum file is empty: {path}")
    return entries


def verify_checksum_file(path, base):
    entries = parse_checksum_file(path)
    for rel, expected in entries.items():
        target = base / rel
        if not target.is_file():
            fail(f"checksum target missing: {target}")
        observed = sha256(target)
        if observed != expected:
            fail(f"checksum mismatch for {target} (observed={observed}, expected={expected})")
    return len(entries)


def compare_checksum_files(observed_path, expected_path):
    observed = parse_checksum_file(observed_path)
    expected = parse_checksum_file(expected_path)
    missing = sorted(set(expected) - set(observed))
    unexpected = sorted(set(observed) - set(expected))
    changed = sorted(k for k in set(observed) & set(expected) if observed[k] != expected[k])
    if missing or unexpected or changed:
        print("FAIL: frozen realistic validation checksum baseline differs", file=sys.stderr)
        for label, items in (("Missing", missing), ("Unexpected", unexpected), ("Changed", changed)):
            if items:
                print(f"  {label} deterministic files:", file=sys.stderr)
                for item in items:
                    print(f"    {item}", file=sys.stderr)
        raise SystemExit(1)
    return len(expected)


if not RUN_DIR.is_dir():
    fail(f"qualification run directory does not exist: {RUN_DIR}")
if not CONFIG_PATH.is_file() or not SAMPLES_PATH.is_file():
    fail(f"run directory must contain copied config.yaml and samples.tsv: {RUN_DIR}")

with CONFIG_PATH.open() as fh:
    RUN_CONFIG = yaml.safe_load(fh)
if not isinstance(RUN_CONFIG, dict):
    fail(f"run config is not a YAML mapping: {CONFIG_PATH}")

# The qualification deliberately teaches the run-owned copy model. Require the copied config to
# use its sibling samples.tsv and output/ directory via explicit absolute paths.
configured_samples = Path(str(RUN_CONFIG.get("samples", ""))).expanduser()
configured_output = Path(str(RUN_CONFIG.get("output", ""))).expanduser()
if not configured_samples.is_absolute():
    fail("qualification config 'samples' must be an absolute path")
if not configured_output.is_absolute():
    fail("qualification config 'output' must be an absolute path")
if configured_samples.resolve() != SAMPLES_PATH:
    fail(f"qualification config must point to its copied samples.tsv: {SAMPLES_PATH}")
expected_output = RUN_DIR / "output"
if configured_output.resolve() != expected_output:
    fail(f"qualification config output must be the run-owned directory: {expected_output}")

if RUN_CONFIG.get("consortium_run") is not True:
    fail("realistic qualification requires consortium_run: true")
ref = RUN_CONFIG.get("reference", {})
expected_reference_identity = {
    "mode": "gdc",
    "genome_fasta": "GRCh38.d1.vd1.fa",
    "gtf": "gencode.v36.annotation.gtf",
    "star_index": "star-2.7.5c_GRCh38.d1.vd1_gencode.v36",
    "sjdb_overhang": 100,
}
for key, expected in expected_reference_identity.items():
    if ref.get(key) != expected:
        fail(f"qualification reference.{key}={ref.get(key)!r}, expected {expected!r}")
reference_dir = Path(str(ref.get("dir", ""))).expanduser()
if not reference_dir.is_absolute():
    fail("qualification reference.dir must be an absolute path to the site GDC installation")
if RUN_CONFIG.get("full_length", {}).get("trim_adapters") is not False:
    fail("realistic qualification requires full_length.trim_adapters: false")
if RUN_CONFIG.get("quantseq", {}).get("bbduk_polyA") is not True:
    fail("realistic qualification requires quantseq.bbduk_polyA: true")
if any(bool(v) for v in RUN_CONFIG.get("modules", {}).values()):
    fail("optional modules must remain disabled for the realistic qualification")
pass_("production-style copied run configuration + consortium GDC identity")

OUTPUT = configured_output.resolve()
RESULTS = OUTPUT / "results"
RESTRICTED = OUTPUT / "restricted"
INTERMEDIATE = OUTPUT / "intermediate"


def validate_pinned_input_semantics():
    manifest = pd.read_csv(TEST / "data_manifest.tsv", sep="\t", dtype=str, keep_default_na=False)
    expected_manifest_files = [
        "SRR493372_1.fastq.gz", "SRR493372_2.fastq.gz", "SRR16932032.fastq.gz"
    ]
    if list(manifest["filename"]) != expected_manifest_files:
        fail("data_manifest.tsv does not contain the intended pinned FASTQ set/order")

    samples = pd.read_csv(SAMPLES_PATH, sep="\t", dtype=str, keep_default_na=False)
    if list(samples["library_id"]) != LIBRARIES:
        fail(f"run sample-sheet library order differs: {list(samples['library_id'])}")
    samples = samples.set_index("library_id")
    fl = samples.loc["REAL_FL"]
    qs = samples.loc["REAL_QS_UMI"]

    if (fl["sample_id"], fl["assay"], fl["layout"], fl["strandedness"], fl["has_umi"], fl["run_accession"]) != (
        "REAL_FL", "full_length", "paired", "unstranded", "false", "SRR493372"
    ):
        fail("REAL_FL sample-sheet metadata differs from the maintained qualification definition")
    if (qs["sample_id"], qs["assay"], qs["layout"], qs["strandedness"], qs["has_umi"], qs["run_accession"]) != (
        "REAL_QS_UMI", "quantseq", "single", "forward", "true", "SRR16932032"
    ):
        fail("REAL_QS_UMI sample-sheet metadata differs from the maintained qualification definition")
    if (qs["umi_pattern"], qs["umi_location"], qs["umi_discard_bases"]) != (
        "NNNNNN", "read1_start", "4"
    ):
        fail("REAL_QS_UMI does not encode the qualified 6-nt UMI + 4-base discard specification")
    if any(fl[col] not in {"", "-"} for col in ("umi_pattern", "umi_location", "umi_discard_bases")):
        fail("REAL_FL must not contain UMI metadata")

    actual_paths = {
        ("REAL_FL", "R1"): Path(fl["fq1"]).expanduser(),
        ("REAL_FL", "R2"): Path(fl["fq2"]).expanduser(),
        ("REAL_QS_UMI", "R1"): Path(qs["fq1"]).expanduser(),
    }
    for key, path in actual_paths.items():
        if not path.is_absolute():
            fail(f"qualification FASTQ path must be absolute for {key}: {path}")
        if not path.is_file():
            fail(f"qualification FASTQ is missing for {key}: {path}")

    # Verify the exact compressed bytes actually referenced by the copied run sheet. This is
    # intentionally repeated here even if get_test_data.sh was used, so a stale/mis-edited path
    # cannot silently pass qualification.
    for _, row in manifest.iterrows():
        key = (row["library_id"], row["read_role"])
        path = actual_paths.get(key)
        if path is None:
            fail(f"sample sheet does not provide the pinned input {key}")
        if path.name != row["filename"]:
            fail(f"{key} references {path.name!r}, expected pinned filename {row['filename']!r}")
        observed_bytes = path.stat().st_size
        if observed_bytes != int(row["bytes"]):
            fail(f"byte-size mismatch for {path}: {observed_bytes} != {row['bytes']}")
        observed_md5 = md5(path)
        if observed_md5 != row["md5"]:
            fail(f"MD5 mismatch for {path}: {observed_md5} != {row['md5']}")
        try:
            with gzip.open(path, "rb") as fh:
                while fh.read(1024 * 1024):
                    pass
        except OSError as exc:
            fail(f"gzip integrity check failed for {path}: {exc}")

    # Explicitly protect the observed Lexogen structure in the exact qualification input.
    # The kit architecture is a 6-nt UMI followed by a 4-nt TATA-like spacer, but the
    # sequenced spacer is not literally TATA in every molecule. The pipeline therefore
    # correctly treats these as four discard bases rather than motif-matched bases.
    #
    # Because the compressed FASTQ identity is already pinned by byte size + MD5 above,
    # this is a structural sanity check rather than a second file-identity test. We scan
    # the full (small) qualification FASTQ so the reported composition is stable and useful
    # when reviewing qualification logs.
    qs_fastq = actual_paths[("REAL_QS_UMI", "R1")]
    expected_spacer = "TATA"
    min_expected_base_fraction = 0.80
    n = 0
    exact_spacer = 0
    lengths = Counter()
    first10 = [Counter() for _ in range(10)]

    with gzip.open(qs_fastq, "rt") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            sequence = fh.readline().rstrip("\n").upper()
            plus = fh.readline()
            quality = fh.readline().rstrip("\n")
            if not plus or not quality:
                fail("SRR16932032 FASTQ is truncated")
            if not header.startswith("@") or not plus.startswith("+"):
                fail("SRR16932032 FASTQ record is malformed")
            if len(sequence) != len(quality):
                fail(f"SRR16932032 qualification read {n + 1} has sequence/quality length mismatch")
            lengths[len(sequence)] += 1
            if len(sequence) < 10:
                fail(f"SRR16932032 qualification read {n + 1} is shorter than the 10-base UMI+discard prefix")
            if sequence[6:10] == expected_spacer:
                exact_spacer += 1
            for pos, base in enumerate(sequence[:10]):
                first10[pos][base] += 1
            n += 1

    if n == 0:
        fail("SRR16932032 FASTQ contains no reads")
    if set(lengths) != {101}:
        fail(f"SRR16932032 qualification reads are not uniformly 101 nt: {dict(lengths)}")

    spacer_fractions = []
    for offset, expected_base in enumerate(expected_spacer, start=6):
        observed = first10[offset][expected_base] / n
        spacer_fractions.append(observed)
        if observed < min_expected_base_fraction:
            fail(
                f"SRR16932032 lacks the expected strong {expected_spacer} spacer signature: "
                f"position {offset + 1} has {expected_base} in only {observed * 100:.2f}% "
                f"of reads (minimum {min_expected_base_fraction * 100:.0f}%)"
            )

    exact_fraction = exact_spacer / n
    umi_composition = []
    for pos in range(6):
        parts = []
        for base in "ACGTN":
            count = first10[pos][base]
            if count:
                parts.append(f"{base}={100.0 * count / n:.2f}%")
        umi_composition.append(f"{pos + 1}:" + ",".join(parts))

    print(
        f"INFO: SRR16932032 reads={n:,}; length=101 nt; exact TATA at positions 7-10="
        f"{exact_spacer:,}/{n:,} ({exact_fraction * 100:.2f}%)"
    )
    print(
        "INFO: SRR16932032 spacer-position support: "
        + ", ".join(
            f"{pos + 7}{base}={fraction * 100:.2f}%"
            for pos, (base, fraction) in enumerate(zip(expected_spacer, spacer_fractions))
        )
    )
    print("INFO: SRR16932032 UMI-position base composition: " + "; ".join(umi_composition))
    pass_("exact pinned ENA bytes + realistic QuantSeq 6-nt UMI / 4-nt TATA-like discard architecture")


validate_pinned_input_semantics()

for directory in (RESULTS, RESTRICTED, INTERMEDIATE):
    if not directory.is_dir():
        fail(f"missing output directory: {directory}")
pass_("results/restricted/intermediate output structure")

# Matrices: exact values will be frozen after qualification; before that, enforce
# deterministic shape/type/non-empty biological signal.
raw_path = RESULTS / "matrices" / "raw_gene_counts.tsv"
raw = pd.read_csv(raw_path, sep="\t").set_index("gene_id")
if list(raw.columns) != LIBRARIES:
    fail(f"raw matrix columns are {list(raw.columns)}, expected {LIBRARIES}")
if len(raw) < 10000:
    fail(f"raw matrix has unexpectedly few genes for GENCODE v36: {len(raw)}")
for col in LIBRARIES:
    values = pd.to_numeric(raw[col], errors="raise")
    if (values < 0).any() or not (values % 1 == 0).all():
        fail(f"raw counts for {col} are not non-negative integers")
    if int(values.sum()) <= 0:
        fail(f"raw counts for {col} contain no assigned expression signal")
raw = raw.astype(int)
pass_("realistic raw GDC gene-count matrix")

umi_path = RESULTS / "matrices" / "umi_molecule_counts.tsv"
umi = pd.read_csv(umi_path, sep="\t").set_index("gene_id")
if list(umi.columns) != UMI_LIBRARIES:
    fail(f"UMI matrix columns are {list(umi.columns)}, expected {UMI_LIBRARIES}")
if list(umi.index) != list(raw.index):
    fail("UMI molecule matrix gene order/set differs from raw matrix")
umi_values = pd.to_numeric(umi["REAL_QS_UMI"], errors="raise")
if (umi_values < 0).any() or not (umi_values % 1 == 0).all() or int(umi_values.sum()) <= 0:
    fail("REAL_QS_UMI molecule matrix is not a non-empty non-negative integer matrix")
umi = umi.astype(int)
pass_("realistic QuantSeq UMI molecule-count matrix")

for library in LIBRARIES:
    path = RESULTS / "libraries" / library / "gene_expression.tsv"
    if not path.is_file():
        fail(f"missing {path}")
    expr = pd.read_csv(path, sep="\t", na_values=["NA"], keep_default_na=True)
    if list(expr.columns) != EXPRESSION_COLUMNS:
        fail(f"{library} gene_expression.tsv has unexpected columns: {list(expr.columns)}")
    expr = expr.set_index("gene_id")
    if list(expr.index) != list(raw.index):
        fail(f"{library} expression gene order/set differs from run raw matrix")
    observed = pd.to_numeric(expr["unstranded"], errors="raise").astype(int)
    if not observed.equals(raw[library]):
        fail(f"{library} expression unstranded counts differ from raw matrix")
    lengths = pd.to_numeric(expr["gene_length"], errors="raise")
    if (lengths <= 0).any():
        fail(f"{library} contains non-positive gene lengths")
    if expr[["gene_name", "gene_type"]].isna().any().any():
        fail(f"{library} has missing GENCODE gene metadata")

    if library == "REAL_FL":
        norm = expr[["fpkm", "fpkm_uq", "tpm"]]
        if norm.isna().any().any():
            fail("REAL_FL normalized expression unexpectedly contains NA")
        if abs(float(pd.to_numeric(expr["tpm"]).sum()) - 1_000_000.0) > 2.0:
            fail("REAL_FL TPM does not sum to approximately 1e6")
        if not expr["umi_molecule_count"].isna().all():
            fail("REAL_FL should have NA umi_molecule_count")
    else:
        if not expr[["fpkm", "fpkm_uq", "tpm"]].isna().all().all():
            fail("REAL_QS_UMI QuantSeq normalized columns should be NA")
        observed_umi = pd.to_numeric(expr["umi_molecule_count"], errors="raise").astype(int)
        if not observed_umi.equals(umi["REAL_QS_UMI"]):
            fail("REAL_QS_UMI per-library molecule counts differ from UMI matrix")
pass_("canonical realistic per-library expression schema and assay-specific semantics")

# Stable QC: avoid arbitrary mapping thresholds, but require coherent non-empty outputs.
qc = pd.read_csv(
    RESULTS / "qc" / "qc_metrics.tsv", sep="\t", na_values=["NA"], keep_default_na=True
).set_index("library_id")
if list(qc.index) != LIBRARIES:
    fail(f"run QC libraries are {list(qc.index)}, expected {LIBRARIES}")
required_umi_qc_columns = {
    "umi_length",
    "umi_location",
    "umi_discard_bases",
    "umi_extract_qc_records",
    "umi_extract_qc_retained_percent",
    "umi_extract_transform_match_percent",
    "umi_extract_tag_match_percent",
}
missing_umi_qc_columns = sorted(required_umi_qc_columns - set(qc.columns))
if missing_umi_qc_columns:
    fail("run QC table is missing UMI extraction QC column(s): " + ", ".join(missing_umi_qc_columns))
expected_meta = {
    "REAL_FL": ("REAL_FL", "full_length", "paired", "false"),
    "REAL_QS_UMI": ("REAL_QS_UMI", "quantseq", "single", "true"),
}
for library in LIBRARIES:
    row = qc.loc[library]
    if (row["sample_id"], row["assay"], row["layout"], str(row["has_umi"]).lower()) != expected_meta[library]:
        fail(f"{library} QC metadata differs from qualification specification")
    for col in ("star_input_records", "uniquely_mapped_reads", "gene_assigned_unstranded"):
        if int(row[col]) <= 0:
            fail(f"{library} QC metric {col} is not positive")
    mapped_pct = float(row["uniquely_mapped_percent"])
    if not 0.0 <= mapped_pct <= 100.0:
        fail(f"{library} uniquely_mapped_percent is outside 0..100")
    per_lib = RESULTS / "libraries" / library / "qc_metrics.tsv"
    if not per_lib.is_file():
        fail(f"missing {per_lib}")
if int(qc.loc["REAL_QS_UMI", "umi_molecules_assigned"]) <= 0:
    fail("REAL_QS_UMI has no assigned UMI molecules")
qs_qc = qc.loc["REAL_QS_UMI"]
if int(qs_qc["umi_length"]) != 6 or str(qs_qc["umi_location"]) != "read1_start" or int(qs_qc["umi_discard_bases"]) != 4:
    fail("REAL_QS_UMI stable QC does not record the qualified 6+4 read1_start UMI design")
if int(qs_qc["umi_extract_qc_records"]) != 10_000:
    fail("REAL_QS_UMI UMI extraction QC did not check the expected deterministic 10,000-record sample")
for col in (
    "umi_extract_qc_retained_percent",
    "umi_extract_transform_match_percent",
    "umi_extract_tag_match_percent",
):
    if float(qs_qc[col]) != 100.0:
        fail(f"REAL_QS_UMI {col} should be 100% for the pinned qualification FASTQ sample")
if not pd.isna(qc.loc["REAL_FL", "umi_molecules_assigned"]):
    fail("REAL_FL non-UMI QC should have NA umi_molecules_assigned")
for col in required_umi_qc_columns:
    if not pd.isna(qc.loc["REAL_FL", col]):
        fail(f"REAL_FL non-UMI QC should have NA {col}")
pass_("realistic per-library and run-level QC metrics + UMI extraction conformance")

# Restricted site-retained products.
for library in LIBRARIES:
    bam = RESTRICTED / "libraries" / library / "alignments" / "genomic.sorted.bam"
    bai = RESTRICTED / "libraries" / library / "alignments" / "genomic.sorted.bam.bai"
    fastqc = RESTRICTED / "libraries" / library / "qc" / "fastqc"
    if not bam.is_file() or bam.stat().st_size == 0 or not bai.is_file():
        fail(f"restricted BAM/index missing for {library}")
    if not fastqc.is_dir() or not list(fastqc.glob("*_fastqc.html")):
        fail(f"restricted FastQC output missing for {library}")
umi_bam = RESTRICTED / "libraries" / "REAL_QS_UMI" / "alignments" / "umi_dedup.bam"
if not umi_bam.is_file() or umi_bam.stat().st_size == 0:
    fail("REAL_QS_UMI deduplicated BAM is missing")
pass_("restricted alignment and detailed QC products")

# Portable run metadata must demonstrate the maintained GDC identity without exposing raw paths.
portable_libraries = pd.read_csv(RESULTS / "run" / "libraries.tsv", sep="\t", dtype=str, keep_default_na=False)
if list(portable_libraries["library_id"]) != LIBRARIES:
    fail("portable libraries.tsv differs from qualification libraries")
if "fq1" in portable_libraries.columns or "fq2" in portable_libraries.columns:
    fail("portable libraries.tsv unexpectedly contains raw FASTQ paths")
qs_row = portable_libraries.set_index("library_id").loc["REAL_QS_UMI"]
if (qs_row["umi_pattern"], qs_row["umi_location"], qs_row["umi_discard_bases"]) != (
    "NNNNNN", "read1_start", "4"
):
    fail("portable REAL_QS_UMI metadata lost the qualified UMI specification")

with (RESULTS / "run" / "config.yaml").open() as fh:
    portable_config = yaml.safe_load(fh)
ref = portable_config.get("reference", {})
if ref != {
    "mode": "gdc",
    "genome_fasta": "GRCh38.d1.vd1.fa",
    "gtf": "gencode.v36.annotation.gtf",
    "star_index": "star-2.7.5c_GRCh38.d1.vd1_gencode.v36",
    "sjdb_overhang": 100,
}:
    fail(f"portable config reference identity differs from maintained GDC bundle: {ref}")
if portable_config.get("full_length", {}).get("trim_adapters") is not False:
    fail("realistic qualification must use GDC-style full_length.trim_adapters=false")
if portable_config.get("quantseq", {}).get("bbduk_polyA") is not True:
    fail("realistic qualification must use quantseq.bbduk_polyA=true")

refs = pd.read_csv(RESULTS / "run" / "references.tsv", sep="\t", dtype=str)
expected_ref_rows = {
    ("genome_fasta", "GRCh38.d1.vd1.fa", "gdc", "NCI GDC GRCh38.d1.vd1"),
    ("annotation_gtf", "gencode.v36.annotation.gtf", "gdc", "NCI GDC GENCODE v36"),
    ("star_index", "star-2.7.5c_GRCh38.d1.vd1_gencode.v36", "gdc", "NCI GDC STAR 2.7.5c pre-built index"),
}
if set(map(tuple, refs[["role", "file", "reference_mode", "source"]].itertuples(index=False, name=None))) != expected_ref_rows:
    fail("references.tsv does not record the intended maintained GDC bundle")
for path in (
    RESTRICTED / "run" / "libraries.original.tsv",
    RESTRICTED / "run" / "config.effective.yaml",
    RESULTS / "run" / "provenance.yaml",
    RESULTS / "run" / "software_versions.tsv",
):
    if not path.is_file():
        fail(f"missing run metadata: {path}")
pass_("portable GDC/run metadata + restricted original run inputs")

software = pd.read_csv(RESULTS / "run" / "software_versions.tsv", sep="\t").set_index("tool_id")
expected_tools = ["bbmap", "fastqc", "htseq", "multiqc", "py", "samtools", "star", "umitools", "snakemake"]
if sorted(software.index.tolist()) != sorted(expected_tools):
    fail(f"software tool set differs: {sorted(software.index.tolist())}")
if (software["version"].astype(str).str.strip() == "").any():
    fail("software_versions.tsv contains an empty version")
if str(software.loc["umitools", "version"]) != "1.1.6":
    fail("realistic qualification requires maintained UMI-tools version 1.1.6")
if str(software.loc["umitools", "container_source"]) != (
    "docker://quay.io/biocontainers/umi_tools:1.1.6--py39hbcbf7aa_0"
):
    fail("realistic qualification UMI-tools container source differs from maintained 1.1.6 image")
pass_("software-version provenance")

multiqc = RESULTS / "qc" / "multiqc_report.html"
if not multiqc.is_file() or multiqc.stat().st_size == 0:
    fail("MultiQC report is missing or empty")
if not (RESTRICTED / "qc" / "multiqc_data").is_dir():
    fail("restricted MultiQC data directory is missing")
html = multiqc.read_text(errors="replace")
rendered = " ".join(unescape(re.sub(r"<[^>]+>", " ", html)).split())
for label in (
    "Library ID", "Sample ID", "STAR input", "Unique mapped %",
    "UMI design", "UMI transform %", "UMI molecules",
):
    if label not in rendered:
        fail(f"MultiQC Library QC Summary label missing: {label}")
for library in LIBRARIES:
    if library not in rendered:
        fail(f"MultiQC report does not contain library identifier {library}")
pass_("portable MultiQC report + site-retained detailed MultiQC data")

# Package integrity and deterministic cross-run/cross-site subset.
package_n = verify_checksum_file(RESULTS / "run" / "checksums.sha256", RESULTS)
validation_path = RESULTS / "run" / "validation_checksums.sha256"
validation_n = verify_checksum_file(validation_path, RESULTS)
pass_(f"package checksums ({package_n} files)")
pass_(f"validation checksums ({validation_n} deterministic files)")

if ARGS.skip_frozen_baseline:
    pass_("frozen realistic cross-installation baseline comparison intentionally skipped")
else:
    frozen = EXPECTED / "validation_checksums.sha256"
    if not frozen.is_file():
        fail(
            "realistic frozen baseline has not been created yet; during initial qualification use "
            "--skip-frozen-baseline, reproduce the generated validation manifest in two clean runs, "
            "then freeze tests/real/expected/validation_checksums.sha256"
        )
    frozen_n = compare_checksum_files(validation_path, frozen)
    if frozen_n != validation_n:
        fail(f"frozen/generated validation entry counts differ ({frozen_n} vs {validation_n})")
    pass_(f"frozen realistic cross-installation validation baseline ({frozen_n} files)")

print("\nRealistic public-data qualification PASSED.")
