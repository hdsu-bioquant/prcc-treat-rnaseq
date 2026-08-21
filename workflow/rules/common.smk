# common.smk — config parsing, library model, path policy, containers, final outputs
import os
import re
import sys
import yaml
from os.path import join

# ---- run output contract ---------------------------------------------------#
# ``output`` is the run root. ``results`` is the portable canonical data product,
# ``restricted`` contains site-retained sequence/infrastructure-sensitive products,
# and ``intermediate`` contains disposable processing artefacts.
OUTPUT_ROOT = config.get("output")
if not OUTPUT_ROOT:
    raise ValueError("Run config must define 'output: /path/to/run/output'")

RESULTS      = join(OUTPUT_ROOT, "results")
RESTRICTED   = join(OUTPUT_ROOT, "restricted")
INTERMEDIATE = join(OUTPUT_ROOT, "intermediate")
TMPDIR       = config.get("tmpdir", join(INTERMEDIATE, "tmp"))

# ---- references ------------------------------------------------------------#
REFDIR    = config["reference"]["dir"]
FASTA     = join(REFDIR, config["reference"]["genome_fasta"])
GTF       = join(REFDIR, config["reference"]["gtf"])
STAR_IDX  = join(REFDIR, config["reference"]["star_index"])
STAR_IDX_DONE = join(STAR_IDX, "SAindex")
GENE_LENGTHS  = join(REFDIR, "gene_lengths.tsv")

# ---- sequencing-library sample sheet --------------------------------------#
# One row = one sequencing library. library_id drives processing/output naming;
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


def library_assay(library):
    return samples.loc[library, "assay"]


UMI_LIBRARIES = [lib for lib in LIBRARIES if library_has_umi(lib)]
FL_UMI_LIBRARIES = [lib for lib in FL_LIBRARIES if library_has_umi(lib)]
QS_UMI_LIBRARIES = [lib for lib in QS_LIBRARIES if library_has_umi(lib)]

LIBRARY_PATTERN = "|".join([re.escape(lib) for lib in LIBRARIES]) if LIBRARIES else "x"
FL_LIBRARY_PATTERN = "|".join([re.escape(lib) for lib in FL_LIBRARIES]) if FL_LIBRARIES else "__no_fl_libraries__"
QS_LIBRARY_PATTERN = "|".join([re.escape(lib) for lib in QS_LIBRARIES]) if QS_LIBRARIES else "__no_qs_libraries__"
FL_UMI_LIBRARY_PATTERN = "|".join([re.escape(lib) for lib in FL_UMI_LIBRARIES]) if FL_UMI_LIBRARIES else "__no_fl_umi_libraries__"
QS_UMI_LIBRARY_PATTERN = "|".join([re.escape(lib) for lib in QS_UMI_LIBRARIES]) if QS_UMI_LIBRARIES else "__no_qs_umi_libraries__"

wildcard_constraints:
    library = LIBRARY_PATTERN


def raw_fq1(wc):
    return samples.loc[wc.library, "fq1"]


def raw_fq2(wc):
    return samples.loc[wc.library, "fq2"]


# ---- generic per-library paths --------------------------------------------#
def _ilib(library, *parts):
    return join(INTERMEDIATE, "libraries", library, *parts)


def _rlib(library, *parts):
    return join(RESTRICTED, "libraries", library, *parts)


def _plib(library, *parts):
    return join(RESULTS, "libraries", library, *parts)


# UMI extraction is independent of assay. QuantSeq then continues into its
# assay-specific trimming, while full-length continues directly (or via fastp).
def qs_trim_input(wc):
    if library_has_umi(wc.library):
        return _ilib(wc.library, "preprocess", "R1.umi.fastq.gz")
    return samples.loc[wc.library, "fq1"]


def fl_pretrim_input(wc):
    if library_has_umi(wc.library):
        return {
            "fq1": _ilib(wc.library, "preprocess", "R1.umi.fastq.gz"),
            "fq2": _ilib(wc.library, "preprocess", "R2.umi.fastq.gz"),
        }
    return {"fq1": samples.loc[wc.library, "fq1"], "fq2": samples.loc[wc.library, "fq2"]}


def fl_star_input(wc):
    if config.get("full_length", {}).get("trim_adapters", False):
        return {
            "fq1": _ilib(wc.library, "preprocess", "R1.trim.fastq.gz"),
            "fq2": _ilib(wc.library, "preprocess", "R2.trim.fastq.gz"),
        }
    return fl_pretrim_input(wc)


HTSEQ_STRAND = {"forward": "yes", "reverse": "reverse", "unstranded": "no"}
def htseq_strand(wc):
    return HTSEQ_STRAND[samples.loc[wc.library, "strandedness"]]


def gene_expression_inputs(wc):
    inp = {
        "counts": _ilib(wc.library, "star", "ReadsPerGene.out.tab"),
        "lengths": GENE_LENGTHS,
    }
    if library_has_umi(wc.library):
        inp["umi_counts"] = _ilib(wc.library, "quantification", "umi_molecule_counts.tsv")
    return inp


# ---- container/software registry ------------------------------------------#
# One workflow-owned manifest is the source of truth for pinned container URIs
# and the software versions reported in provenance / MultiQC. Production runs
# do not add one version-probe job per container; release tests can verify that
# these declarations match the image contents.
SOFTWARE_MANIFEST = os.path.abspath(os.path.join(workflow.basedir, "config", "software_versions.yaml"))
with open(SOFTWARE_MANIFEST) as _fh:
    _software_manifest = yaml.safe_load(_fh)
if not isinstance(_software_manifest, dict) or not isinstance(_software_manifest.get("tools"), dict):
    raise ValueError(f"Invalid software manifest: {SOFTWARE_MANIFEST}")
TOOL_REGISTRY = _software_manifest["tools"]
for _tool_name, _tool_spec in TOOL_REGISTRY.items():
    _missing = {"software", "version", "uri", "pull_default"} - set(_tool_spec)
    if _missing:
        raise ValueError(
            f"Software manifest entry {_tool_name!r} is missing: {', '.join(sorted(_missing))}"
        )

CONTAINER_DIR = os.path.abspath(os.path.join(workflow.basedir, "..", "containers", "sif"))
def _img(name, uri):
    local = os.path.join(CONTAINER_DIR, name + ".sif")
    return local if os.path.exists(local) else uri

IMG = {name: _img(name, spec["uri"]) for name, spec in TOOL_REGISTRY.items()}


def used_tool_names():
    used = {"star", "fastqc", "multiqc", "samtools", "py"}
    if QS_LIBRARIES:
        used.add("bbmap")
    if UMI_LIBRARIES:
        used.update({"umitools", "htseq"})
    if FL_LIBRARIES and config.get("full_length", {}).get("trim_adapters", False):
        used.add("fastp")
    modules = config.get("modules", {})
    if modules.get("rseqc", False):
        used.add("rseqc")
    if modules.get("fusion", False):
        used.update({"arriba", "starfusion"})
    if modules.get("te", False):
        used.add("tetx")
    if modules.get("ase", False):
        used.add("gatk")
    return sorted(used)


# ---- final outputs ---------------------------------------------------------#
def all_outputs(wc):
    out = [FASTA, GTF, STAR_IDX_DONE, GENE_LENGTHS]

    # Canonical portable per-library products.
    out += expand(join(RESULTS, "libraries/{lib}/gene_expression.tsv"), lib=LIBRARIES)
    out += expand(join(RESULTS, "libraries/{lib}/qc_metrics.tsv"), lib=LIBRARIES)

    # Run-level canonical products.
    out += [
        join(RESULTS, "matrices/raw_gene_counts.tsv"),
        join(RESULTS, "qc/qc_metrics.tsv"),
        join(RESULTS, "qc/multiqc_report.html"),
        join(RESULTS, "run/libraries.tsv"),
        join(RESULTS, "run/config.yaml"),
        join(RESULTS, "run/provenance.yaml"),
        join(RESULTS, "run/software_versions.tsv"),
        join(RESULTS, "run/references.tsv"),
        join(RESULTS, "run/manifest.tsv"),
        join(RESULTS, "run/checksums.sha256"),
        join(RESULTS, "run/validation_checksums.sha256"),
    ]
    if UMI_LIBRARIES:
        out += [join(RESULTS, "matrices/umi_molecule_counts.tsv")]

    # Site-retained alignment products. STAR's unsorted and transcriptome BAMs
    # remain disposable intermediates; the coordinate-sorted genomic BAM is kept.
    out += expand(join(RESTRICTED, "libraries/{lib}/alignments/genomic.sorted.bam"), lib=LIBRARIES)
    out += expand(join(RESTRICTED, "libraries/{lib}/alignments/genomic.sorted.bam.bai"), lib=LIBRARIES)
    if UMI_LIBRARIES:
        out += expand(join(RESTRICTED, "libraries/{lib}/alignments/umi_dedup.bam"), lib=UMI_LIBRARIES)

    # Unsanitized run metadata remains site-local by default.
    out += [
        join(RESTRICTED, "run/libraries.original.tsv"),
        join(RESTRICTED, "run/config.effective.yaml"),
    ]

    # Optional modules. Aggregate QC may be portable; variant/fusion-like and
    # other not-yet-reviewed feature outputs remain restricted by default.
    modules = config.get("modules", {})
    if modules.get("rseqc", False):
        out += expand(join(RESULTS, "libraries/{lib}/{lib}.rseqc_read_distribution.txt"), lib=FL_LIBRARIES)
    if modules.get("fusion", False):
        out += expand(join(RESTRICTED, "libraries/{lib}/features/fusion/arriba.fusions.tsv"), lib=FL_LIBRARIES)
        out += expand(join(RESTRICTED, "libraries/{lib}/features/fusion/starfusion.predictions.tsv"), lib=FL_LIBRARIES)
    if modules.get("te", False):
        out += expand(join(RESTRICTED, "libraries/{lib}/features/te/TEcount.cntTable"), lib=FL_LIBRARIES)
    if modules.get("ase", False):
        out += expand(join(RESTRICTED, "libraries/{lib}/features/ase/ASEReadCounter.tsv"), lib=FL_LIBRARIES)

    return out
