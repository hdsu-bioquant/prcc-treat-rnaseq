# pRCC-RNA-Seq

pRCC-RNA-Seq is an assay-aware Snakemake workflow for reproducible bulk RNA-seq processing and gene-expression quantification. The supported core currently covers full-length paired-end and QuantSeq single-end libraries, optional library-level UMI handling, GDC-aligned STAR processing, standardized QC, and portable run outputs.

The pipeline was developed for the **pRCC-TREAT consortium** to support harmonized multi-site processing of papillary renal cell carcinoma RNA-seq data. It can also be used independently outside the consortium, including with deliberate non-consortium reference configurations.

## Workflow overview

```text
        FULL-LENGTH poly-A (PE)                    QUANTSEQ 3' tag (SE)
                raw FASTQ                               raw FASTQ
                    │                                      │
             UMI extract if present                 UMI extract if present
                    │                                      │
             optional fastp                        BBDuk poly(A)/adapter trim
                    │                                      │
                    └──────────────┬───────────────────────┘
                                   ▼
                         STAR two-pass (GDC parameters)
                                   │
                         STAR GeneCounts (raw)
                                   │
                 ┌─────────────────┴─────────────────┐
                 ▼                                   ▼
       canonical gene_expression.tsv       UMI dedup + molecule counts
        (FL also: FPKM/UQ/TPM)                when UMI is present
                 └─────────────────┬─────────────────┘
                                   ▼
                     run matrices + standardized QC
```

The canonical expression basis is STAR unstranded raw `GeneCounts`. Full-length libraries additionally receive GDC-style FPKM, FPKM-UQ, and TPM values. UMI-bearing libraries additionally receive molecule counts after deterministic UMI deduplication.

## Supported core workflow

| Feature | Supported core |
|---|---|
| Full-length RNA-seq | paired-end |
| QuantSeq 3' RNA-seq | single-end |
| Strandedness | unstranded, forward, reverse |
| UMI handling | optional, library-specific |
| Supported UMI architecture | fixed contiguous 5' UMI on R1 or R2, with optional adjacent discard bases |
| Alignment | STAR two-pass with maintained GDC parameters |
| Consortium reference | GRCh38.d1.vd1 + GENCODE v36 + maintained GDC STAR 2.7.5c index |
| Expression outputs | raw counts for all libraries; FPKM/FPKM-UQ/TPM for full-length libraries |
| Execution | local workstation or SLURM HPC |

Optional fusion, TE, ASE, and RSeQC modules are not part of the supported consortium core unless explicitly stated otherwise for a release.

## Requirements and compatibility

The workflow is intended for Linux systems with Conda and Apptainer or Singularity available. Workflow tools such as STAR, samtools, HTSeq, UMI-tools, FastQC, MultiQC, fastp, and BBDuk run inside maintained containers and do not require separate host installation.

| Component | Maintained / practical target |
|---|---|
| Controller environment | `environments/controller.yaml` |
| Python | 3.13.x in the maintained controller environment |
| Snakemake | 9.19.0 |
| Container runtime | Apptainer or Singularity |
| Execution profiles | maintained local and SLURM templates |
| Synthetic qualification | fixed at 2 cores; lightweight installation test |
| Human local execution | must be able to schedule a 64,000 MB STAR job |
| Practical human workstation | 8+ usable cores and approximately 80–96 GB+ RAM recommended |
| Consortium GDC resources | maintained human GDC reference bundle; approximately 59 GB installed |

The 64,000 MB STAR request is a workflow scheduling requirement, not a claim that every complete run has a universal fixed memory minimum. For SLURM execution, appropriate worker nodes or partitions must satisfy the workflow's per-job resource requests.

## Quick start: verify the software installation

The quickest useful first run is the same lightweight first-line installation path used by consortium sites. From the repository root:

```bash
conda env create -f environments/controller.yaml
conda activate prcc-rnaseq-controller

python --version
snakemake --version
apptainer --version        # or: singularity --version

bash containers/pull_images.sh
bash tests/synthetic/run_test.sh
```

If the synthetic test passes, the maintained controller/container stack can execute the supported core workflow against the frozen synthetic validation baseline.

From there:

- **general users:** continue with [`docs/users/general/`](docs/users/general/README.md);
- **consortium sites:** continue with [`SOP-02 — Site qualification`](docs/users/consortium/SOPs/SOP-02-site-qualification.md) to install and qualify the GDC resources, configure a persistent local or SLURM profile, and run the realistic public-data qualification.

For an existing installation, do not recreate the controller environment unnecessarily; use the maintained documentation and qualification procedures appropriate to the task.

## Outputs at a glance

Every run uses a three-part output layout:

```text
output/
├── results/       portable canonical results and sanitized run metadata
├── restricted/    site-retained sequence-level / infrastructure-sensitive products
└── intermediate/  disposable workflow intermediates
```

Canonical portable expression output is written per library as:

```text
results/libraries/<library_id>/gene_expression.tsv
```

Run-level matrices include:

```text
results/matrices/raw_gene_counts.tsv
results/matrices/umi_molecule_counts.tsv   # when UMI libraries are present
```

For consortium runs, `results/` is the intended portable delivery package; `restricted/` remains local unless separately approved. See the [consortium results contract](docs/users/consortium/io-definitions/results-contract.md) or the [general outputs guide](docs/users/general/outputs.md) for details.

## Reproducibility and qualification

The repository maintains explicit controls for reproducible operation:

- pipeline-owned release metadata in `workflow/release.yaml`;
- a maintained controller environment in `environments/controller.yaml`;
- pinned container URIs and expected tool versions in `workflow/config/software_versions.yaml`;
- a canonical GDC installed-reference manifest and site qualification stamp;
- deterministic synthetic qualification;
- realistic public-data qualification against a maintained validation baseline;
- a read-only site preflight in `scripts/verify_installation.py`;
- maintainer-facing static release checks in `tests/release/check_release_consistency.py`.

Exact SIF byte identity across sites is not required. Container presence/version probes plus synthetic and realistic qualification are used to establish operational conformity.

## Documentation

### General users

[`docs/users/general/`](docs/users/general/README.md) provides a conventional pipeline guide covering installation, configuration, execution, outputs, and troubleshooting. General users may deliberately depart from the consortium resource contract using `consortium_run: false` where supported.

### Consortium users

[`docs/users/consortium/`](docs/users/consortium/README.md) contains the controlled, self-contained consortium documentation:

- numbered SOPs under [`docs/users/consortium/SOPs/`](docs/users/consortium/SOPs/);
- input/output definitions and contracts under [`docs/users/consortium/io-definitions/`](docs/users/consortium/io-definitions/).

Consortium sites should follow the applicable SOP/IO document versions for the designated pipeline release.

### Maintainers

Release policy, qualification-baseline maintenance, release checks, and tracked technical debt are documented under [`docs/maintainers/`](docs/maintainers/README.md).

## Repository-maintained operational components

- controller environment: [`environments/controller.yaml`](environments/controller.yaml)
- run templates: [`templates/`](templates/README.md)
- local/SLURM profile templates: [`templates/profiles/`](templates/profiles/README.md)
- software/container manifest: `workflow/config/software_versions.yaml`
- GDC reference resources: [`resources/`](resources/README.md)
- site installation preflight: `scripts/verify_installation.py`
- synthetic/realistic qualification: [`tests/`](tests/README.md)

Pipeline release identity is maintained in `workflow/release.yaml`; copied run configurations do not control it.

## Development and authorship

pRCC-RNA-Seq was developed for the [pRCC-TREAT consortium](https://www.eppermed.eu/funding-projects/projects-results/project-database/prcc-treat/) by the Department of Bioinformatics at the [Institute of Pharmacy and Molecular Biotechnology (IPMB), Heidelberg University](https://www.ipmb.uni-heidelberg.de/en), with development based in the [Health Data Science Unit / BioQuant group](https://www.hdsu.org).

**Software authors**

- **Jan-Eric Bökenkamp** — lead development and maintenance
- **Robert Schwarz** — foundational pipeline development
- **Carl Herrmann** — principal investigator and supervision

Formal software citation metadata is maintained in `CITATION.cff` for releases. License information and release history are maintained in the repository root.
