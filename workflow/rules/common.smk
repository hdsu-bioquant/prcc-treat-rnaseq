# common.smk — config parsing, library model, path policy, containers, final outputs
import csv
import hashlib
import os
import re
import sys
import yaml
from os.path import join

# Resolve workflow-owned paths through the physical filesystem location.  A user
# may launch Snakemake from a logical/symlinked path that is valid on the login
# node but not reproduced inside Apptainer on compute nodes.  Repository-owned
# files must therefore never depend on the shell's logical $PWD.
WORKFLOW_DIR = os.path.realpath(workflow.basedir)
REPO_DIR = os.path.realpath(os.path.join(WORKFLOW_DIR, ".."))
SCRIPT_DIR = os.path.join(WORKFLOW_DIR, "scripts")
CONFIG_DIR = os.path.join(WORKFLOW_DIR, "config")
RESOURCE_DIR = os.path.join(REPO_DIR, "resources")
CONTAINER_DIR = os.path.join(REPO_DIR, "containers", "sif")

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
TMPDIR = config.get("tmpdir", join(INTERMEDIATE, "tmp"))

# Keep per-user runtime/cache state created by containerized tools out of the
# Snakemake working directory. These files are disposable and belong with
# other run-specific temporary data.
RUNTIME_DIR = os.path.abspath(join(TMPDIR, "runtime"))

MPLCONFIGDIR = join(RUNTIME_DIR, "matplotlib")
XDG_CACHE_HOME = join(RUNTIME_DIR, "cache")
JAVA_USER_HOME = join(RUNTIME_DIR, "java-home")

# Set on the host as well as explicitly inside Apptainer containers.
os.environ["MPLCONFIGDIR"] = MPLCONFIGDIR
os.environ["XDG_CACHE_HOME"] = XDG_CACHE_HOME

os.environ["APPTAINERENV_MPLCONFIGDIR"] = MPLCONFIGDIR
os.environ["APPTAINERENV_XDG_CACHE_HOME"] = XDG_CACHE_HOME
os.environ["APPTAINERENV_JAVA_TOOL_OPTIONS"] = f"-Duser.home={JAVA_USER_HOME}"

# ---- consortium/reference policy -------------------------------------------#
# ``consortium_run`` is deliberately explicit. Consortium runs enforce the
# maintained GDC reference identity; non-consortium runs remain free to use the
# pipeline's local-reference mode. Basic GDC structural checks stay enabled for
# every run that selects reference.mode=gdc.
CONSORTIUM_RUN = config.get("consortium_run")
if not isinstance(CONSORTIUM_RUN, bool):
    raise ValueError(
        "Run config must define top-level 'consortium_run: true' or "
        "'consortium_run: false' as a YAML boolean"
    )

REFERENCE_CONFIG = config.get("reference")
if not isinstance(REFERENCE_CONFIG, dict):
    raise ValueError("Run config must define a 'reference:' mapping")

REFERENCE_MODE = REFERENCE_CONFIG.get("mode", "gdc")
if REFERENCE_MODE not in {"gdc", "local"}:
    raise ValueError(
        f"Unsupported reference.mode: {REFERENCE_MODE!r} (expected 'gdc' or 'local')"
    )

for _key in ("dir", "genome_fasta", "gtf", "star_index"):
    if not REFERENCE_CONFIG.get(_key):
        raise ValueError(f"Run config reference.{_key} must be defined")

REFDIR    = REFERENCE_CONFIG["dir"]
FASTA     = join(REFDIR, REFERENCE_CONFIG["genome_fasta"])
GTF       = join(REFDIR, REFERENCE_CONFIG["gtf"])
STAR_IDX  = join(REFDIR, REFERENCE_CONFIG["star_index"])
STAR_IDX_DONE = join(STAR_IDX, "SAindex")

_GDC_RESOURCE_TABLE = os.path.join(RESOURCE_DIR, "gdc_resources.tsv")
_GDC_CANONICAL_MANIFEST = os.path.join(RESOURCE_DIR, "gdc_installed_reference.sha256")
_GDC_QUALIFICATION_STAMP = ".prcc_treat_reference_qualification.tsv"

# Gene lengths are pipeline-derived metadata, not part of the official GDC
# installation. Keep production reference bundles read-only by generating this
# small helper under the run intermediate tree. Local/synthetic mode retains its
# established reference-local path and behavior.
GENE_LENGTHS = (
    join(INTERMEDIATE, "reference/gene_lengths.tsv")
    if REFERENCE_MODE == "gdc"
    else join(REFDIR, "gene_lengths.tsv")
)


def _missing_or_empty_file(path):
    return not os.path.isfile(path) or os.path.getsize(path) == 0


def _validate_preinstalled_gdc_references():
    """Fast structural validation for every pre-installed GDC-mode run."""
    missing = []
    for _label, _path in (("genome FASTA", FASTA), ("annotation GTF", GTF)):
        if _missing_or_empty_file(_path):
            missing.append(f"{_label}: {_path}")

    _star_required = (
        "Genome",
        "SA",
        "SAindex",
        "chrLength.txt",
        "chrName.txt",
        "chrNameLength.txt",
        "chrStart.txt",
        "genomeParameters.txt",
    )
    if not os.path.isdir(STAR_IDX):
        missing.append(f"STAR index directory: {STAR_IDX}")
    else:
        for _name in _star_required:
            _path = join(STAR_IDX, _name)
            if _missing_or_empty_file(_path):
                missing.append(f"STAR index file: {_path}")

    if missing:
        details = "\n".join(f"  - {item}" for item in missing)
        raise ValueError(
            "GDC reference installation is missing or incomplete.\n"
            f"Configured reference.dir: {REFDIR}\n"
            f"Missing required path(s):\n{details}\n\n"
            "Install/repair the reference bundle before analysis. For the maintained "
            "GDC bundle, for example:\n"
            f"  bash resources/get_gdc_references.sh {REFDIR}\n"
            "Normal workflow execution does not download GDC references."
        )


def _expected_consortium_gdc_paths():
    if not os.path.isfile(_GDC_RESOURCE_TABLE):
        raise ValueError(f"Maintained GDC resource table is missing: {_GDC_RESOURCE_TABLE}")

    expected = {}
    with open(_GDC_RESOURCE_TABLE, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            role = row.get("role", "")
            installed = row.get("installed_path", "")
            if role and installed:
                expected[role] = installed

    required = {"genome_fasta", "annotation_gtf", "star_index"}
    missing = required - set(expected)
    if missing:
        raise ValueError(
            f"Maintained GDC resource table lacks role(s): {', '.join(sorted(missing))}"
        )
    return expected


def _sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_qualification_stamp(path):
    values = {}
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t", 1)
            if len(parts) != 2:
                raise ValueError(
                    f"Malformed consortium reference qualification stamp {path!r} "
                    f"at line {lineno}"
                )
            values[parts[0]] = parts[1]
    return values


def _validate_consortium_reference_identity():
    """Enforce consortium-specific reference identity without rehashing the bundle.

    Before the maintainer-owned canonical installed-reference manifest is frozen,
    development consortium runs enforce the exact maintained GDC filenames and
    structural completeness and emit a warning. Once that manifest exists in the
    repository, ``consortium_run: true`` additionally requires a site qualification
    stamp whose manifest identity matches the shipped canonical manifest.
    """
    if not CONSORTIUM_RUN:
        return

    if REFERENCE_MODE != "gdc":
        raise ValueError(
            "consortium_run: true requires reference.mode: gdc. "
            "Use consortium_run: false for local/custom reference analyses."
        )

    expected = _expected_consortium_gdc_paths()
    configured = {
        "genome_fasta": REFERENCE_CONFIG["genome_fasta"],
        "annotation_gtf": REFERENCE_CONFIG["gtf"],
        "star_index": REFERENCE_CONFIG["star_index"],
    }
    mismatches = [
        f"reference.{('gtf' if role == 'annotation_gtf' else role)}: "
        f"configured {configured[role]!r}, expected {expected[role]!r}"
        for role in ("genome_fasta", "annotation_gtf", "star_index")
        if configured[role] != expected[role]
    ]

    try:
        sjdb_overhang = int(REFERENCE_CONFIG.get("sjdb_overhang", 100))
    except (TypeError, ValueError):
        sjdb_overhang = None
    if sjdb_overhang != 100:
        mismatches.append(
            f"reference.sjdb_overhang: configured {REFERENCE_CONFIG.get('sjdb_overhang')!r}, "
            "expected 100"
        )

    if mismatches:
        details = "\n".join(f"  - {item}" for item in mismatches)
        raise ValueError(
            "consortium_run: true requires the maintained pRCC-TREAT GDC reference identity.\n"
            f"Mismatch(es):\n{details}\n"
            "Set consortium_run: false for a deliberately non-consortium/custom-reference run."
        )

    if not os.path.isfile(_GDC_CANONICAL_MANIFEST):
        print(
            "WARNING: consortium_run=true, but the maintainer-owned canonical installed-reference "
            "SHA256 manifest has not yet been frozen. This is expected during development; "
            "exact GDC filenames and fast structural checks are enforced, but byte-level site "
            "qualification is not yet required.",
            file=sys.stderr,
        )
        return

    if _missing_or_empty_file(_GDC_CANONICAL_MANIFEST):
        raise ValueError(
            f"Canonical consortium reference manifest is empty: {_GDC_CANONICAL_MANIFEST}"
        )

    canonical_sha256 = _sha256_file(_GDC_CANONICAL_MANIFEST)
    stamp_path = join(REFDIR, _GDC_QUALIFICATION_STAMP)
    if _missing_or_empty_file(stamp_path):
        raise ValueError(
            "consortium_run: true requires this installed reference bundle to be qualified "
            "against the maintainer-frozen SHA256 manifest.\n"
            f"Missing qualification stamp: {stamp_path}\n"
            "Run once after installation/copying:\n"
            f"  bash resources/verify_gdc_references.sh --qualify {REFDIR}"
        )

    stamp = _read_qualification_stamp(stamp_path)
    observed = stamp.get("canonical_manifest_sha256")
    if observed != canonical_sha256:
        raise ValueError(
            "The installed GDC reference qualification stamp does not match the canonical "
            "manifest shipped with this pipeline checkout.\n"
            f"Installed stamp manifest SHA256: {observed or '<missing>'}\n"
            f"Pipeline canonical manifest SHA256: {canonical_sha256}\n"
            "Re-qualify the installed bundle with:\n"
            f"  bash resources/verify_gdc_references.sh --qualify {REFDIR}"
        )


if REFERENCE_MODE == "gdc":
    _validate_preinstalled_gdc_references()
_validate_consortium_reference_identity()

# ---- sequencing-library sample sheet --------------------------------------#
# One row = one sequencing library. library_id drives processing/output naming;
# sample_id identifies the underlying biological sample and may repeat.
# Import workflow-owned helpers from the same canonical script directory used by
# containerized rule shell commands below.
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)
from sample_sheet import load_and_validate_samples
from umi import compile_umitools_extract, parse_umi_spec

samples = load_and_validate_samples(config["samples"])

LIBRARIES = list(samples["library_id"])
FL_LIBRARIES = list(samples[samples["assay"] == "full_length"]["library_id"])
QS_LIBRARIES = list(samples[samples["assay"] == "quantseq"]["library_id"])


def library_has_umi(library):
    return samples.loc[library, "has_umi"] == "true"


def biological_sample_id(library):
    return samples.loc[library, "sample_id"]


def library_assay(library):
    return samples.loc[library, "assay"]


def library_layout(library):
    return samples.loc[library, "layout"]


UMI_LIBRARIES = [lib for lib in LIBRARIES if library_has_umi(lib)]
PAIRED_UMI_LIBRARIES = [lib for lib in UMI_LIBRARIES if library_layout(lib) == "paired"]
SINGLE_UMI_LIBRARIES = [lib for lib in UMI_LIBRARIES if library_layout(lib) == "single"]

# Sample-sheet validation has already guaranteed these fields are valid. Compile them
# once here so workflow rules consume pipeline-native UMI specifications rather than
# reinterpreting raw sample-sheet strings independently.
UMI_SPECS = {
    lib: parse_umi_spec(
        samples.loc[lib, "umi_pattern"],
        samples.loc[lib, "umi_location"],
        samples.loc[lib, "umi_discard_bases"],
    )
    for lib in UMI_LIBRARIES
}
UMITOOLS_EXTRACT_SPECS = {lib: compile_umitools_extract(spec) for lib, spec in UMI_SPECS.items()}


def library_umi_spec(library):
    return UMI_SPECS[library]


def library_umitools_extract_spec(library):
    return UMITOOLS_EXTRACT_SPECS[library]


LIBRARY_PATTERN = "|".join([re.escape(lib) for lib in LIBRARIES]) if LIBRARIES else "x"
FL_LIBRARY_PATTERN = "|".join([re.escape(lib) for lib in FL_LIBRARIES]) if FL_LIBRARIES else "__no_fl_libraries__"
QS_LIBRARY_PATTERN = "|".join([re.escape(lib) for lib in QS_LIBRARIES]) if QS_LIBRARIES else "__no_qs_libraries__"
UMI_LIBRARY_PATTERN = "|".join([re.escape(lib) for lib in UMI_LIBRARIES]) if UMI_LIBRARIES else "__no_umi_libraries__"
PAIRED_UMI_LIBRARY_PATTERN = "|".join([re.escape(lib) for lib in PAIRED_UMI_LIBRARIES]) if PAIRED_UMI_LIBRARIES else "__no_paired_umi_libraries__"
SINGLE_UMI_LIBRARY_PATTERN = "|".join([re.escape(lib) for lib in SINGLE_UMI_LIBRARIES]) if SINGLE_UMI_LIBRARIES else "__no_single_umi_libraries__"

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

def umi_primary_input(wc):
    spec = library_umitools_extract_spec(wc.library)
    return samples.loc[wc.library, "fq1" if spec.umi_read == 1 else "fq2"]


def umi_mate_input(wc):
    spec = library_umitools_extract_spec(wc.library)
    return samples.loc[wc.library, "fq2" if spec.umi_read == 1 else "fq1"]


def umi_primary_output(wc):
    spec = library_umitools_extract_spec(wc.library)
    read = "R1" if spec.umi_read == 1 else "R2"
    return _ilib(wc.library, "preprocess", f"{read}.umi.fastq.gz")


def umi_mate_output(wc):
    spec = library_umitools_extract_spec(wc.library)
    read = "R2" if spec.umi_read == 1 else "R1"
    return _ilib(wc.library, "preprocess", f"{read}.umi.fastq.gz")


def umi_dedup_paired_flag(wc):
    return "--paired" if library_layout(wc.library) == "paired" else ""


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
SOFTWARE_MANIFEST = os.path.join(CONFIG_DIR, "software_versions.yaml")
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

def _img(name, spec):
    # A manifest entry may override the historical <tool>.sif filename. The
    # pull helper version-probes tools that declare version_probe before an
    # existing local image is reused.
    local_name = spec.get("local_sif", name + ".sif")
    local = os.path.join(CONTAINER_DIR, local_name)
    return local if os.path.exists(local) else spec["uri"]

IMG = {name: _img(name, spec) for name, spec in TOOL_REGISTRY.items()}


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
