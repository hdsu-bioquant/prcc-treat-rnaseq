# pRCC-RNA-Seq run templates

This directory defines the maintained executable run interface. General guidance is under `docs/users/general/`; the controlled consortium contracts are under `docs/users/consortium/`.

Start every new production analysis by copying `config.yaml` and `samples.tsv` to a
run-specific location. Do not edit the files in `templates/` in place.

The workflow intentionally has no repository-wide biological run configuration. A run is
fully specified by its own copied configuration and library sheet, passed explicitly to
Snakemake with `--configfile`.

## 1. Start a new run

From the repository root:

```bash
mkdir -p /path/to/run
cp templates/config.yaml  /path/to/run/config.yaml
cp templates/samples.tsv  /path/to/run/samples.tsv
```

Edit both copied files, then validate the library sheet before launching anything:

```bash
python workflow/scripts/sample_sheet.py /path/to/run/samples.tsv
```

A successful validation prints the number of sequencing libraries and biological samples.
Then perform a dry-run:

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --configfile /path/to/run/config.yaml \
  -n
```

For a real execution, use one of the official execution-profile templates in
[`templates/profiles/`](profiles/README.md). Copy the appropriate profile once for the local
compute environment, edit the copy, and reuse it across runs. For example:

```bash
mkdir -p ~/.config/snakemake
cp -r templates/profiles/local ~/.config/snakemake/prcc-rnaseq-local
# or: cp -r templates/profiles/slurm ~/.config/snakemake/prcc-rnaseq-slurm

snakemake \
  --snakefile workflow/Snakefile \
  --configfile /path/to/run/config.yaml \
  --profile ~/.config/snakemake/prcc-rnaseq-local \
  --keep-going
```

Execution profiles contain machine/site settings such as local CPU/RAM limits, SLURM
submission limits, partition/account names, latency handling, and optional Apptainer bind
paths. They are deliberately separate from this biological run interface and are not copied
into the portable results package.

## 2. Path rules

**Important:** relative paths are resolved relative to the process working directory
(normally the pRCC-RNA-Seq repository root), **not** relative to `config.yaml` or
`samples.tsv`.

For production runs, absolute paths are recommended for:

- the copied `samples.tsv` path in `config.yaml`;
- the run `output` directory;
- FASTQ paths in `samples.tsv`;
- optional module resources;
- `reference.dir` if GDC resources are stored outside the repository.

The default template uses `reference.dir: resources/gdc`, which is convenient when the
pipeline is launched from the repository root and the standard reference bundle is stored
there.

## 3. `samples.tsv`: one row = one sequencing library

The columns that are always required are:

```text
library_id  sample_id  assay  layout  strandedness  fq1  fq2  has_umi
```

If at least one row has `has_umi=true`, the sheet must also contain all three UMI-detail
columns:

```text
umi_pattern  umi_location  umi_discard_bases
```

If every row has `has_umi=false`, those three columns may be omitted completely. The official
template includes them so it can be used for mixed runs. The two template rows are
**examples only**; replace or delete them, and do not treat their UMI values as protocol
defaults.

| Column | Requirement | Meaning / allowed values |
|---|---|---|
| `library_id` | required, unique | Sequencing-library identifier and workflow/output key. |
| `sample_id` | required | Biological sample identifier. May repeat when multiple libraries belong to the same biological sample. |
| `assay` | required | `full_length` or `quantseq`. |
| `layout` | required | `paired` or `single`. Currently implemented: `full_length + paired`, `quantseq + single`. |
| `strandedness` | required | `unstranded`, `forward`, or `reverse`; set from the library-preparation protocol. |
| `fq1` | required | R1 or single-end FASTQ: `.fastq`, `.fastq.gz`, `.fq`, `.fq.gz`. |
| `fq2` | conditional | Required for paired libraries. Use `-` (or empty) for single-end libraries. |
| `has_umi` | required | `true` or `false`. UMI handling is library-specific and independent of assay. |
| `umi_pattern` | conditional | Pipeline-native fixed-length UMI description. Currently one or more `N` characters only, e.g. `NNNNNN` or `NNNNNNNN`. Required when `has_umi=true`. |
| `umi_location` | conditional | `read1_start` or `read2_start`. `read2_start` requires a paired-end library. Required when `has_umi=true`. |
| `umi_discard_bases` | conditional | Non-negative integer number of additional bases immediately after the extracted UMI that must be removed but are not part of the molecular identifier. Required when `has_umi=true`; `0` means none. |

### Identifier rules

`library_id` and `sample_id` currently enforce only technical filesystem/wildcard safety:

```text
^[A-Za-z0-9][A-Za-z0-9._-]*$
```

That means letters, numbers, `.`, `_`, and `-` are accepted, and the identifier must start
with a letter or number. This is **not** the final consortium semantic naming convention.
Patient/sample/library naming standards can be tightened later without changing the workflow
model.

### Additional metadata columns

Extra columns such as `batch`, `site`, `patient_id`, or `model_id` are accepted and preserved
in the loaded run table. They do not currently change workflow routing.

They are also **not copied automatically into the portable `results/run/libraries.tsv`**.
The portable manifest intentionally contains only the standardized technical library fields;
the exact original sheet remains under `restricted/run/libraries.original.tsv`.

### Validation behavior

The sample sheet is validated while the DAG is parsed. The workflow fails early for, among
other things:

- duplicate or unsafe `library_id` values;
- missing/unsafe `sample_id` values;
- unsupported assay/layout/strandedness values;
- paired libraries without R2 or single-end libraries with an R2 path;
- missing FASTQs or unsupported FASTQ suffixes;
- missing UMI-detail columns when UMI-bearing libraries are present;
- incomplete, contradictory, or unsupported UMI metadata;
- `read2_start` on a single-end library;
- assigning the same FASTQ to more than one library row, including symlink-equivalent paths.

## 4. `config.yaml`

### Maintainer-owned release identity

Pipeline release identity is not a run configuration field. It is maintained in
`workflow/release.yaml` and is injected into portable run metadata by the workflow. Users therefore
cannot accidentally relabel a copied run configuration as another pipeline release.

### Fields users normally edit

`consortium_run` must be an explicit YAML boolean. Keep `consortium_run: true` for pRCC-TREAT
consortium runs. Set it to `false` only for deliberate non-consortium/custom-reference use. In
consortium mode the workflow enforces the maintained GDC reference filenames,
`reference.sjdb_overhang: 100`, and a site qualification stamp matching the frozen canonical
installed-reference SHA256 manifest. No ~30 GB reference hash is recomputed at normal run startup.

| Field | Required? | Guidance |
|---|---|---|
| `consortium_run` | yes | `true` for pRCC-TREAT consortium runs; `false` for deliberate non-consortium/custom-reference analyses. |
| `samples` | yes | Absolute path to the copied run `samples.tsv` is recommended. |
| `output` | yes | Run-wise output root. Must be different for independent runs. |
| `tmpdir` | no | Optional STAR scratch location. If omitted: `<output>/intermediate/tmp`. |
| `reference.dir` | normally | Location of the exact GDC resource bundle. `resources/gdc` is the repository-root default. |
| `star.threads` | site-dependent | Execution/resource setting; may be tuned without changing scientific parameters. |

### Harmonization settings users normally should not change

For consortium-harmonized production runs, retain the template values for:

- `reference.mode: gdc`;
- GDC filenames / STAR-index directory name;
- `reference.sjdb_overhang: 100`;
- `star.gdc_params`;
- `full_length.trim_adapters: false` unless a deliberately non-standard analysis has been agreed;
- `quantseq.bbduk_polyA: true`.

`reference.mode: local` exists for development/test reference bundles and is not the normal
production setting.

### GDC references

The standard production reference is:

- genome: `GRCh38.d1.vd1.fa`;
- annotation: `gencode.v36.annotation.gtf`;
- STAR index: `star-2.7.5c_GRCh38.d1.vd1_gencode.v36`.

These files must already exist under `reference.dir` before a production run starts. Snakemake
never downloads the GDC production bundle. Install and MD5-verify the maintained bundle once at
each site with:

```bash
bash resources/get_gdc_references.sh resources/gdc
```

The official archive MD5s verify downloaded archives, while the frozen canonical
`resources/gdc_installed_reference.sha256` verifies the qualified extracted installation. Consortium
sites qualify an installed or copied bundle once with
`bash resources/verify_gdc_references.sh --qualify /path/to/gdc`; normal runs check the resulting
small stamp rather than hashing the full reference installation repeatedly.

Pinned URLs and official GDC archive MD5s are installation metadata in
`resources/gdc_resources.tsv`; they are intentionally absent from copied run configs. A quick
site check is available with:

```bash
bash resources/verify_gdc_references.sh resources/gdc
```

For multi-site harmonization, use the fixed pre-built resources rather than rebuilding the human
STAR index locally. If references live on a site-specific filesystem root, keep the scientific
reference identity in the run config and put any required Apptainer bind root in the copied
execution profile.

## 5. Assay and UMI behavior

The supported production combinations are currently:

```text
full_length + paired
quantseq    + single
```

UMI presence is independent of assay. A library with `has_umi=true` is routed through UMI
extraction before its assay-specific preprocessing, then through UMI deduplication after
alignment. The current supported UMI class is deliberately narrow: one fixed-length,
contiguous UMI at the 5′ start of R1 or R2, optionally followed by a fixed number of
additional bases to discard. For example:

```text
NNNNNN    read1_start    4
NNNNNNNN  read1_start    0
NNNNNNNN  read2_start    8
```

`umi_discard_bases` describes only how many adjacent bases to remove; it does not encode or
validate their sequence. Partners provide these semantic fields rather than UMI-tools regexes.
Other architectures (internal/3′/split/variable-length UMIs, for example) are not currently
implemented and should not be represented by stretching these fields beyond their defined
meaning.

The canonical expression basis for both assays is STAR **unstranded raw GeneCounts**.
`gene_expression.tsv` also retains both stranded STAR diagnostic count columns. Full-length
libraries receive FPKM, FPKM-UQ and TPM; these fields are `NA` for QuantSeq. UMI-bearing
libraries additionally receive `umi_molecule_count`.

## 6. Optional modules

All optional modules are off in the template because they are outside the core expression
harmonization contract.

| Config | Scope | Additional resource when enabled | Default output class |
|---|---|---|---|
| `modules.rseqc` | full-length | none beyond standard annotation conversion | `results/` QC |
| `modules.fusion` | full-length | `ctat_genome_lib` for STAR-Fusion | `restricted/` |
| `modules.te` | full-length | `te_gtf` | `restricted/` |
| `modules.ase` | full-length | `ase_germline_vcf` | `restricted/` |

Do not enable a module without also supplying its required resource path. Fusion, TE and ASE
outputs remain site-retained by default pending explicit data-sharing decisions.

## 7. Output contract

Every normal run creates:

```text
<output>/
├── results/        # canonical portable result package
├── restricted/     # retained sequence-level / infrastructure-sensitive products
└── intermediate/   # disposable processing artefacts
```

The main portable products are:

```text
results/
├── libraries/<library_id>/
│   ├── gene_expression.tsv
│   └── qc_metrics.tsv
├── matrices/
│   ├── raw_gene_counts.tsv
│   └── umi_molecule_counts.tsv       # only when UMI libraries exist
├── qc/
│   ├── qc_metrics.tsv
│   └── multiqc_report.html
└── run/
    ├── libraries.tsv
    ├── config.yaml
    ├── provenance.yaml
    ├── software_versions.tsv
    ├── references.tsv
    ├── manifest.tsv
    ├── checksums.sha256
    └── validation_checksums.sha256
```

For pRCC-TREAT partner runs, `results/` is the intended central delivery package. The generic
name is deliberate: outside the consortium it simply represents the portable canonical run
result.

For UMI-bearing libraries, the stable `qc_metrics.tsv` rows also report the declared UMI
length/location/discard specification plus a deterministic sampled extraction-conformance
check (sampled retention, exact raw-to-extracted transformation, and read-name UMI tagging).
This confirms that the configured UMI specification was applied; it is not a protocol-specific
claim about UMI randomness or spacer sequence identity.

`checksums.sha256` verifies transfer/package integrity. `validation_checksums.sha256` contains
the deliberately deterministic subset used for cross-installation harmonization testing.
The synthetic installation test compares that subset to the frozen reference baseline.

## 8. What to copy for each new run

A clean run directory can be as small as:

```text
RUN_DIRECTORY/
├── config.yaml
└── samples.tsv
```

The large FASTQs and GDC resources may live elsewhere and be referenced by absolute path.
The workflow writes all analysis products beneath the configured `output` root.

Before production use at a new site, first run:

```bash
bash tests/synthetic/run_test.sh
```

and confirm that the frozen synthetic validation baseline passes. Then follow the realistic
public-data/GDC-resource qualification under `tests/real/`. That exercise deliberately uses the
same operational model described here: configure a site execution profile once, create a fresh
run directory, copy maintained run files into it, edit only run/site paths, and invoke Snakemake
directly with the persistent profile.
