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
GENE_LENGTHS  = join(REFDIR, "gencode.v36.gene_lengths.tsv")     # only for OPTIONAL FPKM/TPM
# 3' exonic-window GTF (+ BED) for the Branch-A 3'-restricted secondary.
# GTF = HTSeq counting features; BED = samtools pre-filter regions (so HTSeq only
# processes the few % of reads that fall in the 3' windows -> fast on full-depth BAMs).
THREEP_GTF    = join(REFDIR, "gencode.v36.3prime_%dbp.gtf" % config["cross_assay"]["window_bp"])
THREEP_BED    = join(REFDIR, "gencode.v36.3prime_%dbp.bed" % config["cross_assay"]["window_bp"])

# ---- sample sheet ----------------------------------------------------------#
samples = pd.read_csv(config["samples"], sep="\t", dtype=str).set_index("sample", drop=False)
samples = samples.fillna("-")
SAMPLES    = list(samples["sample"])
FL_SAMPLES = list(samples[samples["assay"] == "full_length_pe"]["sample"])
QS_SAMPLES = list(samples[samples["assay"] == "quantseq_3prime_se"]["sample"])

wildcard_constraints:
    sample = "|".join([re.escape(s) for s in SAMPLES]) if SAMPLES else "x"

def raw_fq1(wc):  return samples.loc[wc.sample, "fq1"]
def raw_fq2(wc):  return samples.loc[wc.sample, "fq2"]

# QuantSeq read routing (pRCC-TREAT scheme).
# PRIMARY: raw R1 -> (if has_umi) umi_tools extract -> BBDuk -> STAR SE -> STAR GeneCounts (NO dedup).
# SECONDARY (has_umi only): sort -> umi_tools dedup -> HTSeq-count -> UMI-dedup matrix.
def qs_trim_input(wc):
    if config["quantseq"]["has_umi"]:
        return join(RESULTS, "quantseq", wc.sample, wc.sample + ".umi.fastq.gz")
    return samples.loc[wc.sample, "fq1"]

# HTSeq -s value for the UMI-dedup SECONDARY. QuantSeq FWD is truly "forward" (yes);
# the STAR PRIMARY matrix stays unstranded (scheme choice) but keeps the stranded columns.
QS_HTSEQ_STRAND = {"forward": "yes", "reverse": "reverse", "unstranded": "no"}[config["quantseq"]["strandedness"]]

# Full-length STAR input: trimmed (fastp) if enabled, else raw FASTQs (GDC-exact)
def fl_star_input(wc):
    if config["full_length"]["trim_adapters"]:
        return {"fq1": join(RESULTS, "full_length", wc.sample, "trim", wc.sample + "_R1.trim.fastq.gz"),
                "fq2": join(RESULTS, "full_length", wc.sample, "trim", wc.sample + "_R2.trim.fastq.gz")}
    return {"fq1": samples.loc[wc.sample, "fq1"], "fq2": samples.loc[wc.sample, "fq2"]}

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
    # DE: an image providing DESeq2 + edgeR + limma + apeglm (build from containers/de.Dockerfile)
    "de":        _img("de",        "docker://bioconductor/bioconductor_docker:RELEASE_3_18"),
}

# assays with at least one sample, restricted to those configured for DE
DE_ASSAYS = [a for a in config["de"]["assays"]
             if a in set(samples["assay"])] if config["de"]["enabled"] else []

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
    # SECONDARY: QuantSeq UMI-deduplicated HTSeq matrix (only when UMIs present)
    if config["quantseq"]["has_umi"] and QS_SAMPLES:
        out += expand(join(RESULTS, "quantseq/{s}/{s}.dedup_htseq_counts.tsv"), s=QS_SAMPLES)
        out += [join(RESULTS, "matrix/quantseq_umidedup_matrix.tsv")]
    # OPTIONAL GDC-style normalized outputs (FPKM/FPKM-UQ/TPM) — off = matches scheme
    if config["full_length"].get("compute_fpkm_tpm", False):
        out += expand(join(RESULTS, "full_length/{s}/{s}.augmented_star_gene_counts.tsv"), s=FL_SAMPLES)
    # SECONDARY (Branch A): 3'-restricted HTSeq gene-count matrix
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
