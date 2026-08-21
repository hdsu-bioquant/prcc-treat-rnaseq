# common.smk — config parsing, sample model, container registry, final outputs
import os
import re
import pandas as pd
from os.path import join

# ---- paths -----------------------------------------------------------------#
RESULTS   = config["results"]
TMPDIR    = config["tmpdir"]
REFDIR    = config["reference"]["dir"]
FASTA     = join(REFDIR, config["reference"]["genome_fasta"])
GTF       = join(REFDIR, config["reference"]["gtf"])
STAR_IDX  = join(REFDIR, config["reference"]["star_index"])     # directory
STAR_IDX_DONE = join(STAR_IDX, "SAindex")                        # build/extract marker
GENE_LENGTHS  = join(REFDIR, "gene_lengths.tsv")     # only for OPTIONAL FPKM/TPM

# Legacy 3'-window paths. Kept temporarily so the current workflow still parses
# while the obsolete cross-assay experiment is removed in a later cleanup.
CROSS_ASSAY = config.get("cross_assay", {})
THREEP_WINDOW_BP = int(CROSS_ASSAY.get("window_bp", 500))
THREEP_GTF = join(REFDIR, "gencode.v36.3prime_%dbp.gtf" % THREEP_WINDOW_BP)
THREEP_BED = join(REFDIR, "gencode.v36.3prime_%dbp.bed" % THREEP_WINDOW_BP)

# ---- sample sheet ----------------------------------------------------------#
samples = pd.read_csv(config["samples"], sep="\t", dtype=str).fillna("-")

# UMI metadata now belongs to each library/sample. For backwards compatibility
# with existing consortium sheets, missing columns inherit the old QuantSeq-wide
# settings; new sheets should provide has_umi/umi_pattern/umi_location explicitly.
qcfg = config.get("quantseq", {})
if "has_umi" not in samples.columns:
    samples["has_umi"] = "false"
    if qcfg.get("has_umi", False):
        samples.loc[samples["assay"] == "quantseq_3prime_se", "has_umi"] = "true"
if "umi_pattern" not in samples.columns:
    samples["umi_pattern"] = "-"
    old_pattern = str(qcfg.get("umi_pattern", "-"))
    umi_qs = (samples["assay"] == "quantseq_3prime_se") & samples["has_umi"].astype(str).str.lower().isin(["true", "1", "yes", "y"])
    samples.loc[umi_qs, "umi_pattern"] = old_pattern
if "umi_location" not in samples.columns:
    samples["umi_location"] = "-"
    has_any_umi = samples["has_umi"].astype(str).str.lower().isin(["true", "1", "yes", "y"])
    samples.loc[has_any_umi, "umi_location"] = "read1_start"

samples = samples.set_index("sample", drop=False)
SAMPLES    = list(samples["sample"])
FL_SAMPLES = list(samples[samples["assay"] == "full_length_pe"]["sample"])
QS_SAMPLES = list(samples[samples["assay"] == "quantseq_3prime_se"]["sample"])

_TRUE = {"true", "1", "yes", "y"}
def sample_has_umi(sample):
    return str(samples.loc[sample, "has_umi"]).strip().lower() in _TRUE

def sample_umi_pattern(sample):
    return str(samples.loc[sample, "umi_pattern"]).strip()

def sample_umi_location(sample):
    return str(samples.loc[sample, "umi_location"]).strip()

UMI_SAMPLES    = [s for s in SAMPLES if sample_has_umi(s)]
FL_UMI_SAMPLES = [s for s in FL_SAMPLES if sample_has_umi(s)]
QS_UMI_SAMPLES = [s for s in QS_SAMPLES if sample_has_umi(s)]

# This first generic UMI implementation supports the consortium's current UMI
# structure: a fixed-length UMI at the start of R1. Other locations can be added
# explicitly later rather than silently interpreting them incorrectly.
for s in UMI_SAMPLES:
    if sample_umi_pattern(s) in {"", "-"}:
        raise ValueError("Sample %s has has_umi=true but no umi_pattern" % s)
    if sample_umi_location(s) != "read1_start":
        raise ValueError("Sample %s uses unsupported umi_location=%s (currently supported: read1_start)" %
                         (s, sample_umi_location(s)))

wildcard_constraints:
    sample = "|".join([re.escape(s) for s in SAMPLES]) if SAMPLES else "x"

def raw_fq1(wc):  return samples.loc[wc.sample, "fq1"]
def raw_fq2(wc):  return samples.loc[wc.sample, "fq2"]

# UMI extraction is independent of assay. QuantSeq then continues into its
# assay-specific poly(A)/adapter trimming, while full-length continues directly
# (or through fastp if that option is enabled).
def qs_trim_input(wc):
    if sample_has_umi(wc.sample):
        return join(RESULTS, "quantseq", wc.sample, wc.sample + ".umi.fastq.gz")
    return samples.loc[wc.sample, "fq1"]

def fl_pretrim_input(wc):
    if sample_has_umi(wc.sample):
        return {
            "fq1": join(RESULTS, "full_length", wc.sample, "umi", wc.sample + "_R1.umi.fastq.gz"),
            "fq2": join(RESULTS, "full_length", wc.sample, "umi", wc.sample + "_R2.umi.fastq.gz"),
        }
    return {"fq1": samples.loc[wc.sample, "fq1"], "fq2": samples.loc[wc.sample, "fq2"]}

# Full-length STAR input: trimmed (fastp) if enabled, else UMI-extracted/raw FASTQs.
def fl_star_input(wc):
    if config["full_length"]["trim_adapters"]:
        return {"fq1": join(RESULTS, "full_length", wc.sample, "trim", wc.sample + "_R1.trim.fastq.gz"),
                "fq2": join(RESULTS, "full_length", wc.sample, "trim", wc.sample + "_R2.trim.fastq.gz")}
    return fl_pretrim_input(wc)

HTSEQ_STRAND = {"forward": "yes", "reverse": "reverse", "unstranded": "no"}
def htseq_strand(wc):
    strand = str(samples.loc[wc.sample, "strandedness"]).strip().lower()
    if strand not in HTSEQ_STRAND:
        raise ValueError("Unsupported strandedness for %s: %s" % (wc.sample, strand))
    return HTSEQ_STRAND[strand]

# ---- container registry ---------------------------------------------------#
# Prefer a PRE-PULLED local image at containers/sif/<name>.sif (run
# containers/pull_images.sh once) and fall back to the docker:// URI otherwise.
# Pre-pulling decouples the run from flaky registry access (e.g. quay.io TLS timeouts).
CONTAINER_DIR = os.path.abspath(os.path.join(workflow.basedir, "..", "containers", "sif"))
def _img(name, uri):
    local = os.path.join(CONTAINER_DIR, name + ".sif")
    return local if os.path.exists(local) else uri
IMG = {
    "star":      _img("star",      "docker://quay.io/biocontainers/star:2.7.5c--0"),
    "fastqc":    _img("fastqc",    "docker://quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"),
    "multiqc":   _img("multiqc",   "docker://quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0"),
    "fastp":     _img("fastp",     "docker://quay.io/biocontainers/fastp:0.23.4--hadf994f_2"),
    "bbmap":     _img("bbmap",     "docker://quay.io/biocontainers/bbmap:39.06--h92535d8_0"),
    "umitools":  _img("umitools",  "docker://quay.io/biocontainers/umi_tools:1.1.4--py39hf95cd2a_2"),
    "samtools":  _img("samtools",  "docker://quay.io/biocontainers/samtools:1.19--h50ea8bc_0"),
    "subread":   _img("subread",   "docker://quay.io/biocontainers/subread:2.0.6--he4a0461_2"),
    "htseq":     _img("htseq",     "docker://quay.io/biocontainers/htseq:2.0.9--py39h918f1d6_0"),
    "rseqc":     _img("rseqc",     "docker://quay.io/biocontainers/rseqc:5.0.3--py39hf95cd2a_0"),
    "arriba":    _img("arriba",    "docker://quay.io/biocontainers/arriba:1.1.0--h2e03b76_3"),
    "starfusion":_img("starfusion","docker://trinityctat/starfusion:1.6.0"),
    "tetx":      _img("tetx",      "docker://quay.io/biocontainers/tetranscripts:2.2.3--pyhdfd78af_0"),
    "gatk":      _img("gatk",      "docker://broadinstitute/gatk:4.5.0.0"),
    "py":        _img("py",        "docker://quay.io/biocontainers/pandas:1.5.2"),
    # DE entry kept temporarily until the dedicated DE cleanup pass.
    "de":        _img("de",        "docker://bioconductor/bioconductor_docker:RELEASE_3_18"),
}

# Legacy DE bookkeeping kept parse-safe until the dedicated cleanup pass.
decfg = config.get("de", {})
DE_ASSAYS = [a for a in decfg.get("assays", []) if a in set(samples["assay"])] if decfg.get("enabled", False) else []

# ---- final outputs ---------------------------------------------------------#
def all_outputs(wc):
    out = []
    # references
    out += [FASTA, GTF, STAR_IDX_DONE]
    # PRIMARY (both branches): STAR raw gene counts (unstranded + stranded diagnostics)
    out += expand(join(RESULTS, "full_length/{s}/{s}.star_gene_counts.tsv"), s=FL_SAMPLES)
    out += expand(join(RESULTS, "quantseq/{s}/{s}.star_gene_counts.tsv"), s=QS_SAMPLES)
    # PRIMARY: cohort raw-count matrix (STAR unstranded, both branches) + unified QC
    out += [join(RESULTS, "matrix/gene_counts_matrix.tsv"),
            join(RESULTS, "qc/multiqc_report.html")]
    # SECONDARY: UMI-deduplicated molecule-level counts, independent of assay.
    if UMI_SAMPLES:
        out += expand(join(RESULTS, "full_length/{s}/{s}.dedup_htseq_counts.tsv"), s=FL_UMI_SAMPLES)
        out += expand(join(RESULTS, "quantseq/{s}/{s}.dedup_htseq_counts.tsv"), s=QS_UMI_SAMPLES)
        out += [join(RESULTS, "matrix/umi_dedup_matrix.tsv")]
    # OPTIONAL GDC-style normalized outputs (FPKM/FPKM-UQ/TPM)
    if config["full_length"].get("compute_fpkm_tpm", False):
        out += expand(join(RESULTS, "full_length/{s}/{s}.augmented_star_gene_counts.tsv"), s=FL_SAMPLES)
    # optional modules (full-length only)
    if config["modules"]["rseqc"]:
        out += expand(join(RESULTS, "full_length/{s}/qc/{s}.rseqc.read_distribution.txt"), s=FL_SAMPLES)
    if config["modules"]["fusion"]:
        out += expand(join(RESULTS, "full_length/{s}/fusion/{s}.arriba.fusions.tsv"), s=FL_SAMPLES)
        out += expand(join(RESULTS, "full_length/{s}/fusion/{s}.starfusion.predictions.tsv"), s=FL_SAMPLES)
    if config["modules"]["te"]:
        out += expand(join(RESULTS, "full_length/{s}/te/{s}.TEcount.cntTable"), s=FL_SAMPLES)
    if config["modules"]["ase"]:
        out += expand(join(RESULTS, "full_length/{s}/ase/{s}.ASEReadCounter.tsv"), s=FL_SAMPLES)
    return out
