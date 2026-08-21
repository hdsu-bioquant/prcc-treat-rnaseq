# common.smk — config parsing, library model, container registry, final outputs
import os
import re
import sys
from os.path import join

# ---- paths -----------------------------------------------------------------#
RESULTS   = config["results"]
TMPDIR    = config["tmpdir"]
REFDIR    = config["reference"]["dir"]
FASTA     = join(REFDIR, config["reference"]["genome_fasta"])
GTF       = join(REFDIR, config["reference"]["gtf"])
STAR_IDX  = join(REFDIR, config["reference"]["star_index"])     # directory
STAR_IDX_DONE = join(STAR_IDX, "SAindex")                        # build/extract marker
GENE_LENGTHS  = join(REFDIR, "gene_lengths.tsv")                  # only for OPTIONAL FPKM/TPM

# ---- sequencing-library sample sheet --------------------------------------#
# One row = one sequencing library. library_id drives workflow/output naming;
# sample_id identifies the underlying biological sample and may repeat.
_SCRIPT_DIR = os.path.abspath(os.path.join(workflow.basedir, "scripts"))
if _SCRIPT_DIR not in sys.path:
    sys.path.insert(0, _SCRIPT_DIR)
from sample_sheet import load_and_validate_samples

samples = load_and_validate_samples(config["samples"])

LIBRARIES = list(samples["library_id"])
FL_LIBRARIES = list(samples[samples["assay"] == "full_length"]["library_id"])
QS_LIBRARIES = list(samples[samples["assay"] == "quantseq"]["library_id"])


def library_has_umi(library):
    return samples.loc[library, "has_umi"] == "true"


def library_umi_pattern(library):
    return samples.loc[library, "umi_pattern"]


def library_umi_location(library):
    return samples.loc[library, "umi_location"]


def biological_sample_id(library):
    return samples.loc[library, "sample_id"]


UMI_LIBRARIES = [lib for lib in LIBRARIES if library_has_umi(lib)]
FL_UMI_LIBRARIES = [lib for lib in FL_LIBRARIES if library_has_umi(lib)]
QS_UMI_LIBRARIES = [lib for lib in QS_LIBRARIES if library_has_umi(lib)]

wildcard_constraints:
    library = "|".join([re.escape(lib) for lib in LIBRARIES]) if LIBRARIES else "x"


def raw_fq1(wc):
    return samples.loc[wc.library, "fq1"]


def raw_fq2(wc):
    return samples.loc[wc.library, "fq2"]


# UMI extraction is independent of assay. QuantSeq then continues into its
# assay-specific poly(A)/adapter trimming, while full-length continues directly
# (or through fastp if that option is enabled).
def qs_trim_input(wc):
    if library_has_umi(wc.library):
        return join(RESULTS, "quantseq", wc.library, wc.library + ".umi.fastq.gz")
    return samples.loc[wc.library, "fq1"]


def fl_pretrim_input(wc):
    if library_has_umi(wc.library):
        return {
            "fq1": join(RESULTS, "full_length", wc.library, "umi", wc.library + "_R1.umi.fastq.gz"),
            "fq2": join(RESULTS, "full_length", wc.library, "umi", wc.library + "_R2.umi.fastq.gz"),
        }
    return {"fq1": samples.loc[wc.library, "fq1"], "fq2": samples.loc[wc.library, "fq2"]}


# Full-length STAR input: trimmed (fastp) if enabled, else UMI-extracted/raw FASTQs.
def fl_star_input(wc):
    if config["full_length"]["trim_adapters"]:
        return {
            "fq1": join(RESULTS, "full_length", wc.library, "trim", wc.library + "_R1.trim.fastq.gz"),
            "fq2": join(RESULTS, "full_length", wc.library, "trim", wc.library + "_R2.trim.fastq.gz"),
        }
    return fl_pretrim_input(wc)


HTSEQ_STRAND = {"forward": "yes", "reverse": "reverse", "unstranded": "no"}
def htseq_strand(wc):
    return HTSEQ_STRAND[samples.loc[wc.library, "strandedness"]]


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
}


# ---- final outputs ---------------------------------------------------------#
def all_outputs(wc):
    out = []
    # references
    out += [FASTA, GTF, STAR_IDX_DONE]
    # PRIMARY (both branches): STAR raw gene counts (unstranded + stranded diagnostics)
    out += expand(join(RESULTS, "full_length/{lib}/{lib}.star_gene_counts.tsv"), lib=FL_LIBRARIES)
    out += expand(join(RESULTS, "quantseq/{lib}/{lib}.star_gene_counts.tsv"), lib=QS_LIBRARIES)
    # PRIMARY: cohort raw-count matrix (STAR unstranded, both branches) + unified QC
    out += [join(RESULTS, "matrix/gene_counts_matrix.tsv"),
            join(RESULTS, "qc/multiqc_report.html")]
    # SECONDARY: UMI-deduplicated molecule-level counts, independent of assay.
    if UMI_LIBRARIES:
        out += expand(join(RESULTS, "full_length/{lib}/{lib}.dedup_htseq_counts.tsv"), lib=FL_UMI_LIBRARIES)
        out += expand(join(RESULTS, "quantseq/{lib}/{lib}.dedup_htseq_counts.tsv"), lib=QS_UMI_LIBRARIES)
        out += [join(RESULTS, "matrix/umi_dedup_matrix.tsv")]
    # OPTIONAL GDC-style normalized outputs (FPKM/FPKM-UQ/TPM)
    if config["full_length"].get("compute_fpkm_tpm", False):
        out += expand(join(RESULTS, "full_length/{lib}/{lib}.augmented_star_gene_counts.tsv"), lib=FL_LIBRARIES)
    # optional modules (full-length only)
    modules = config.get("modules", {})
    if modules.get("rseqc", False):
        out += expand(join(RESULTS, "full_length/{lib}/qc/{lib}.rseqc.read_distribution.txt"), lib=FL_LIBRARIES)
    if modules.get("fusion", False):
        out += expand(join(RESULTS, "full_length/{lib}/fusion/{lib}.arriba.fusions.tsv"), lib=FL_LIBRARIES)
        out += expand(join(RESULTS, "full_length/{lib}/fusion/{lib}.starfusion.predictions.tsv"), lib=FL_LIBRARIES)
    if modules.get("te", False):
        out += expand(join(RESULTS, "full_length/{lib}/te/{lib}.TEcount.cntTable"), lib=FL_LIBRARIES)
    if modules.get("ase", False):
        out += expand(join(RESULTS, "full_length/{lib}/ase/{lib}.ASEReadCounter.tsv"), lib=FL_LIBRARIES)
    return out
