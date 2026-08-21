
# RNA-Seq-pRCC pipeline — User guide

HDSU University Heidelberg 2026
RNA-Seq pipeline for the pRCC-TREAT consortium
Prof. Dr. Carl Herrmann, Dr. Jan-Eric Bökenkamp, Robert Schwarz, B.Sc.

## Executive Summary

- unified, **GDC-aligned, assay-aware** RNA-seq pipeline for the multi-site papillary
renal-cell carcinoma (pRCC) harmonization effort (**pRCC-TREAT**)
- standardized to the NCI GDC mRNA Analysis pipeline (GRCh38.d1.vd1 + GENCODE v36 + STAR 2.7.5c + STAR GeneCounts)
- normalisation + cross-protocol integration are downstream
- primary output: raw-count matrix

```
   FULL-LENGTH poly-A (PE)                        QUANTSEQ 3' tag (SE + UMI)
   raw FASTQ                                      raw FASTQ
     │  (optional) fastp                            │  umi_tools extract  (6 bp UMI)
     │                                              │  BBDuk  (polyA + adapter right-trim)
     ▼                                              ▼
   STAR 2-pass (GDC params) ◀─ GDC STAR index ─▶ STAR 2-pass (SE, GDC params)
     │                                              │
     ▼  STAR --quantMode GeneCounts                 ▼  STAR --quantMode GeneCounts (non-dedup)
   raw gene counts (unstranded = PRIMARY)         raw gene counts (unstranded = PRIMARY)
     │  (optional) FPKM / FPKM-UQ / TPM             │  umi_tools dedup ─▶ HTSeq  (SECONDARY)
     └───────────────────────┬──────────────────────┘
                             ▼
             harmonized raw-count matrix  +  unified MultiQC
```
> The pipeline is on GitHub: **https://github.com/hdsu-bioquant/pRCC-RNA-Seq**.
> The pipeline runs anywhere that provides **Conda + Snakemake + Apptainer/Singularity** — an HPC cluster with a scheduler (e.g., SLURM) or a single workstation.
> It does not depend on any particular site, filesystem, or authentication system.

## Table of contents

0. Initial Setup
1. Session Setup
2. Sample Sheet construction
3. Pipeline configuration
4. Pipeline execution
5. Outputs & Interpretation
6. Optional modules
7. Reproducibility — containers and provenance
8. References

---

## 0. Initial Setup

You only need this section the first time you ever use the pipeline. Subsequent sessions skip to §1.

### 0.1 Requirements

| Software | Explanatory Section |
|---|---|
| **Conda / Miniconda** | §0.2 |
| **Snakemake 8 or 9** | §0.3 |
| **Apptainer or Singularity** | §0.4 |
| **GDC reference files** (genome + GTF + STAR index) | §0.5 |
| **Container images** (pre-pulled) | §0.6 |
| **Sample Sheet + FASTQ files** | §2 |

The following resources are fetched automatically by the pipeline (large):

| Resource | Size | Notes |
|---|---|---|
| **resources/gdc/** | ~59 GB | GRCh38.d1.vd1 genome + GENCODE v36 GTF + GDC-built STAR 2.7.5c index. Downloaded once (§0.5). |
| **containers/sif/** | ~1.4 GB | Pinned biocontainers pulled once to local `.sif` (§0.6). |

**The following software does *not* need to be installed:** Docker (Apptainer reads Docker images),
STAR, samtools, HTSeq, UMI-tools, fastp, BBDuk, FastQC, MultiQC — these are all inside containers. The
reference genome/annotation/index are downloaded from the GDC (§0.5).


### 0.2 Miniconda Installation (skip if already available)

```bash
curl -fsSL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh
bash /tmp/miniconda.sh -b -p $HOME/miniconda3
rm -f /tmp/miniconda.sh
$HOME/miniconda3/bin/conda init bash
source ~/.bashrc
conda --version
```

### 0.3 Create the Snakemake pipeline environment

```bash
conda create -y -n snakemake-9.19.0 \
    -c conda-forge -c bioconda \
    'snakemake=9.19.0' \
    snakemake-executor-plugin-slurm
```

Notes:
- `snakemake-executor-plugin-slurm` is **required** to dispatch jobs via `--profile profiles/slurm`. Without it, the SLURM profile cannot submit jobs. On a different scheduler, install the matching plugin instead (e.g. `snakemake-executor-plugin-cluster-generic`, `-lsf`, `-drmaa`).
- Snakemake 8 is supported as well (`'snakemake>=8'`).

Verify:
```bash
conda activate snakemake-9.19.0
snakemake --version
```

### 0.4 Make Apptainer / Singularity available

Every rule runs inside a container via Apptainer (or Singularity). You do **not** need to install it if your
cluster already provides it. Load whichever name your site uses:

```bash
module avail                       # find the apptainer / singularity module on your cluster
module load apptainer              # or whatever your site calls it, e.g.: module load singularity
apptainer --version                # or: singularity --version
```

If your cluster does not use environment modules, just make sure `apptainer` (or `singularity`) is on your
`PATH`. On a personal workstation, install Apptainer once from your OS package manager.

### 0.5 Get the GDC reference files (one-time, ~59 GB)

Two equivalent options — both land the files in `resources/gdc/`:

```bash
# (a) let Snakemake fetch them automatically (config: reference.download_references: true — the default), OR
# (b) pre-download + MD5-verify yourself:
bash resources/get_gdc_references.sh resources/gdc
```

This retrieves the **exact GDC files** (MD5-checked): `GRCh38.d1.vd1.fa`,
`gencode.v36.annotation.gtf`, and the GDC-built `star-2.7.5c_GRCh38.d1.vd1_gencode.v36` index.
The index is **downloaded, never rebuilt**, so it is byte-identical to GDC's / TCGA's.

### 0.6 Pre-pull the container images (one-time, ~1.4 GB)

```bash
# Apptainer/Singularity must be available first (see §0.4)
bash containers/pull_images.sh
```

Pre-pulling to `containers/sif/<name>.sif` decouples runs from flaky registry access (e.g. quay.io TLS
timeouts). The pipeline uses the local `.sif` if present, else falls back to the `docker://` URI (pulled
on demand into your Apptainer cache).

---

## 1. Session Setup

### 1.1 Connect and start a tmux session

An SSH disconnect kills every process whose terminal closes — including a multi-hour controller. Run inside
`tmux` so the Snakemake controller survives disconnects.

```bash
ssh your_username@your-cluster-login-node

tmux new -s rnaseq     # or `tmux attach -t rnaseq` to re-attach
```

### 1.2 Make Apptainer / Singularity available

```bash
module load apptainer     # your cluster's module name (see §0.4); skip if already on PATH
```

### 1.3 Snakemake environment activation

```bash
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate snakemake-9.19.0
```

### 1.4 Navigation to the pipeline directory

```bash
cd /path/to/pRCC-RNA-Seq
```

> **Warning:** the directory path must NOT contain `&` or whitespace.

---

## 2. Sample Sheet construction

The sample sheet is a **tab-separated file with a header row**. One row represents one
**sequencing library**. The run-specific config points to this file via `samples:`.

Required columns:

```
library_id  sample_id  assay  layout  strandedness  fq1  fq2  has_umi  umi_pattern  umi_location
```

| Column | Meaning / allowed values |
|---|---|
| `library_id` | Unique sequencing-library identifier. This is the workflow/output identifier and must be unique. |
| `sample_id` | Biological sample identifier. It may occur in multiple rows when one biological sample has multiple libraries. |
| `assay` | `full_length` \| `quantseq` |
| `layout` | `paired` \| `single`. Currently implemented combinations are `full_length + paired` and `quantseq + single`. |
| `strandedness` | `unstranded` \| `forward` \| `reverse` |
| `fq1` | Path to R1/single-end FASTQ (`.fastq`, `.fastq.gz`, `.fq`, `.fq.gz`). |
| `fq2` | Path to R2 for paired-end libraries; use `-` or leave empty for single-end libraries. |
| `has_umi` | `true` \| `false`. UMI handling is a **library property**, independent of assay. |
| `umi_pattern` | UMI-tools extraction pattern, e.g. `NNNNNN`; required when `has_umi=true`, otherwise `-`. |
| `umi_location` | Currently `read1_start` for UMI libraries; otherwise `-`. |

Additional metadata columns (for example `batch`, `site`, `patient_id`, `model_id`) are
allowed and preserved in the loaded table, but do not currently control workflow routing.
Formal consortium naming conventions for patient/sample/library identifiers will be defined
separately; the pipeline currently enforces only filesystem-safe identifiers (letters,
numbers, `.`, `_`, `-`).

Example:

```
library_id  sample_id  assay        layout  strandedness  fq1                    fq2                    has_umi  umi_pattern  umi_location  batch
FL_01       S01        full_length  paired  reverse       /abs/S01_R1.fastq.gz   /abs/S01_R2.fastq.gz   false    -            -             run1
QS_01       S01        quantseq     single  forward       /abs/S01_QS.fastq.gz   -                      true     NNNNNN       read1_start   run2
```

The sheet is validated before the DAG is built. The pipeline fails immediately for, among
other things, duplicate `library_id` values, missing FASTQs, invalid assay/layout/strand
values, paired libraries without R2, inconsistent UMI fields, unsupported UMI locations,
or accidental reuse of the same FASTQ in multiple library rows.

---


## 3. Pipeline configuration

Each analysis run supplies its own YAML configuration explicitly with `--configfile`.
There is intentionally **no repository-wide default run configuration**: omitted run-specific
settings should fail rather than silently inherit values from an unrelated example dataset.
A documented `templates/` directory will be added once the user-facing schema is frozen.

### 3.1 Main settings

```yaml
samples: /path/to/run/samples.tsv       # library sheet (§2)
results: /path/to/run/results           # output destination
tmpdir:  /path/to/run/results/tmp       # STAR scratch

reference:
  dir: resources/gdc
  genome_fasta: GRCh38.d1.vd1.fa
  gtf:          gencode.v36.annotation.gtf
  star_index:   star-2.7.5c_GRCh38.d1.vd1_gencode.v36
  sjdb_overhang: 100
```

### 3.2 Assay-level settings

UMI presence, UMI structure, read layout, and strandedness are **not** QuantSeq-global
configuration switches. They are specified per library in the sample sheet.

```yaml
full_length:
  trim_adapters:    false
  count_column:     unstranded
  compute_fpkm_tpm: true

quantseq:
  bbduk_polyA: true
```

The primary output is STAR raw gene counts on the non-deduplicated read basis for both
assays. If a library has `has_umi=true`, UMI extraction happens before assay-specific
processing and a UMI-deduplicated molecule-count layer is additionally produced.

| | Full-length | QuantSeq 3′ |
|---|---|---|
| Current supported layout | paired-end | single-end |
| UMI | optional, library-level | optional, library-level |
| Trim | optional fastp; off by default | BBDuk poly(A)/adapter trim |
| Align | STAR two-pass GDC parameters | STAR two-pass GDC parameters |
| Primary counts | STAR GeneCounts, non-deduplicated | STAR GeneCounts, non-deduplicated |
| Secondary UMI counts | when `has_umi=true` | when `has_umi=true` |
| Optional normalization | FPKM/FPKM-UQ/TPM | none in pipeline |

The primary cohort matrix uses the STAR **unstranded** column for both assays for uniformity
with TCGA. Per-library STAR tables retain the stranded diagnostic columns, while the
sample-sheet `strandedness` controls HTSeq counting of UMI-deduplicated BAMs.

### 3.3 Tunable STAR parameters

```yaml
star:
  threads: 8
  gdc_params: >-
    --twopassMode Basic --quantMode TranscriptomeSAM GeneCounts
    --outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNoverLmax 0.1
    --alignSJoverhangMin 8 --chimSegmentMin 15  [... full GDC block ...]
```

For production GDC adherence, the canonical GDC STAR recipe should not be changed.

## 4. Pipeline execution

Ideal start status: inside `tmux`, Apptainer/Singularity available, `snakemake-9.19.0` active, working
directory at the pipeline root, and references (§0.5) + images (§0.6) present.

### 4.1 Dry-run

```bash
snakemake -n --configfile /path/to/run/config.yaml
```

`-n` parses the sample sheet, builds the DAG, and prints what would run — without scheduling jobs or
pulling containers. Read the job-count summary that is printed at the bottom of the output after the command has run.

### 4.2 Launching the pipeline

Depending on the **infrastructure** of your working environment, the **snakemake workflow** can be edited to use a **profile** personalised to your architecture: **HPC cluster** or **local working station** (see the steps below):

**On an HPC cluster (SLURM):**

```bash
snakemake --profile profiles/slurm \
          --configfile /path/to/run/config.yaml \
          --rerun-incomplete --keep-going
```
- `--profile profiles/slurm` reads `profiles/slurm/config.yaml` (SLURM executor, apptainer-only deployment, `jobs: 50`, `rerun-triggers: mtime`).
- `--rerun-incomplete` redoes any rule whose outputs have been marked as incomplete.
- `--keep-going` — when one sample's job fails, keep scheduling the others.

> **Adapt the SLURM profile to your cluster** — edit `profiles/slurm/config.yaml`:
> - **partition/queue:** `default-resources` sets `slurm_partition=single`. Change `single` to your cluster's partition name (and add `slurm_account=<your_account>` if your site requires an account).
> - **container bind path:** `apptainer-args: "--bind <path>"` must point at a directory that contains your FASTQs, the references, and the results, so containers can read/write them. Set it to the parent of your data (e.g. `--bind /scratch/<you>` or several `--bind` paths).
> - optionally tune `jobs`, `latency-wait`, and `default-resources` (`mem_mb`, `runtime`) for your queue limits.
>
> A non-SLURM scheduler (SGE/LSF/PBS) needs the matching `snakemake-executor-plugin-*` (§0.3) and an equivalent profile.


**On a local workstation:**

```bash
snakemake --profile profiles/local \
          --configfile /path/to/run/config.yaml \
          --rerun-incomplete --keep-going
```
`profiles/local` uses `cores: 16, jobs: 16` — tune to your machine, and set its `apptainer-args` bind to your
paths. The human STAR index needs **~30 GB RAM**, so keep concurrent STAR jobs low (e.g. one at a time on a
96 GB box).

> **Long runs:** keep the controller inside `tmux` (§1.1), or submit the Snakemake controller itself as a
> scheduler job so it survives login-node reboots or SSH disconnects.

### 4.3 Monitoring progress

```bash
squeue -u "$USER"                                             # SLURM queue (adjust for your scheduler)
tail -f .snakemake/log/$(ls -t .snakemake/log/*.snakemake.log | head -1)   # controller log
ls logs/slurm/                                               # per-rule SLURM logs
```

### 4.4 Run locally without a scheduler

> This is good for small tests or debugging.

```bash
snakemake --cores 4 --use-singularity --configfile /path/to/run/config.yaml
```

Skip `--profile profiles/slurm`. Only for small tests — STAR alignment of the human genome needs ~30 GB RAM.

---

---

## 5. Outputs & Interpretation

Under `results/`:

```
qc/multiqc_report.html                                aggregated QC (FastQC + STAR)
qc/fastqc/<library>.done                               per-library FastQC marker (+ html/zip)

full_length/<s>/
  ├── <s>.Aligned.sortedByCoord.bam(.bai)             coord-sorted genome BAM
  ├── <s>.ReadsPerGene.out.tab                        raw STAR GeneCounts (3 strand columns)
  ├── <s>.star_gene_counts.tsv                        ★ PRIMARY per-library counts (cleaned)
  └── <s>.augmented_star_gene_counts.tsv              FPKM / FPKM-UQ / TPM   (if compute_fpkm_tpm)

quantseq/<s>/
  ├── <s>.star_gene_counts.tsv                        ★ PRIMARY QuantSeq counts (non-dedup, uniform basis)
  ├── <s>.dedup.bam                                   UMI-deduplicated BAM (secondary path)
  └── <s>.dedup_htseq_counts.tsv                      SECONDARY UMI-dedup HTSeq counts

matrix/
  ├── gene_counts_matrix.tsv                          ★★ PRIMARY cohort matrix (STAR unstranded, BOTH branches)
  └── umi_dedup_matrix.tsv                         SECONDARY UMI-dedup matrix (all UMI libraries)
```
`matrix/gene_counts_matrix.tsv` is the **primary output** — raw counts, gene × library, with `# sample_id`, `# assay`, and `# layout` annotation rows.

### 5.1 Library QC guide

| Metric | Source | Expected healthy values |
|---|---|---|
| Uniquely mapped % | `*/Log.final.out` (in MultiQC) | full-length ≳ 90 %, QuantSeq ≳ 85 % |
| % reads assigned to genes | STAR `N_*` rows / MultiQC | protocol-dependent |
| TPM sum per library | `augmented_star_gene_counts.tsv` | ≈ 1e6 (sanity check) |
| Duplication rate | FastQC (in MultiQC) | high for QuantSeq is **expected** (3'-tag) |
| Library size (assigned counts) | matrix column sums | full-length typically 3–4× deeper than QuantSeq |

A library failing on mapping % or gene-assignment fraction usually has an upstream problem (degraded RNA,
adapter/contamination, wrong reference) that no downstream step can fix.

---

## 6. Optional modules

### 6.1 FPKM / FPKM-UQ / TPM (GDC augmented output)

```yaml
full_length:
  compute_fpkm_tpm: true    
```
Output:
`full_length/<s>/<s>.augmented_star_gene_counts.tsv`.

This module reproduces GDC's `augmented_star_gene_counts`: `FPKM = RCg·1e9/(RCpc·L)`, `FPKM-UQ` (75th-percentile
protein-coding denominator), `TPM = (RCg/L)/Σ(RCj/L)·1e6` and is **only appliable to full-length data**. 
3′-tag data has no length bias, so FPKM/TPM are invalid for QuantSeq (CPM needs to be used downstream).

### 6.2 UMI-deduplicated secondary

UMI handling is configured per library in the sample sheet, for example:

```
library_id  sample_id  assay        layout  strandedness  fq1       fq2       has_umi  umi_pattern  umi_location
LIB01       S01        full_length  paired  reverse       R1.fq.gz  R2.fq.gz  true     NNNNNN       read1_start
LIB02       S02        quantseq     single  forward       QS.fq.gz  -         true     NNNNNN       read1_start
```

For `has_umi=true`, `umi_tools extract` moves the UMI into the read name before alignment.
After alignment, `umi_tools dedup` collapses PCR duplicates by mapping position + UMI and
HTSeq produces molecule-level gene counts. The secondary cohort output is
`matrix/umi_dedup_matrix.tsv` and may contain UMI-bearing libraries from either assay.

---

## 7. Reproducibility — containers and provenance

Every rule declares a `container:` image, and the profiles run **Apptainer** only:

```yaml
software-deployment-method:
  - apptainer
apptainer-args: "--bind /path/that/contains/your/data/refs/results"   # set this for your site (see §4.2)
```

| Pinned image | Used by |
|---|---|
| `quay.io/biocontainers/star:2.7.5c--0` | STAR alignment + GeneCounts (both branches) |
| `quay.io/biocontainers/samtools:1.19--h50ea8bc_0` | sort / index |
| `quay.io/biocontainers/umi_tools:1.1.4--py39hf95cd2a_2` | UMI extract + dedup (UMI-bearing libraries) |
| `quay.io/biocontainers/bbmap:39.06--h92535d8_0` | BBDuk polyA/adapter trim (QuantSeq) |
| `quay.io/biocontainers/fastp:0.23.4--hadf994f_2` | optional full-length adapter trimming |
| `quay.io/biocontainers/htseq:2.0.9--py39h918f1d6_0` | HTSeq counting (UMI-dedup secondary) |
| `quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0` · `multiqc:1.21--pyhdfd78af_0` | QC + aggregation |
| `quay.io/biocontainers/pandas:1.5.2` | count extraction / matrix merge / FPKM-TPM |

Reference: **GRCh38.d1.vd1 + GENCODE v36 + GDC STAR 2.7.5c index**, all GDC-exact and MD5-verified.
`_img()` prefers a local `containers/sif/<name>.sif` (pull once with `containers/pull_images.sh`, §0.6) and
falls back to the `docker://` URI otherwise.

The profiles set `rerun-triggers: mtime` — Snakemake re-runs a job only when its **input files** change (by
mtime), not when a rule's code, parameters, or container reference changes. Editing comments or config
therefore does not trigger a full rerun.

---

## 8. References

### Workflow infrastructure
- **Snakemake.** Mölder F, *et al.* F1000Research 2021;10:33.
- **Apptainer / Singularity.** Kurtzer GM *et al.* PLoS ONE 2017;12:e0177459.

### GDC standard + reference
- **GDC mRNA Analysis pipeline.** NCI Genomic Data Commons (DR32+); docs.gdc.cancer.gov · github.com/NCI-GDC/gdc-rnaseq-cwl.
- **GENCODE.** Frankish A *et al.* NAR 2021;49:D916.

### Alignment
- **STAR.** Dobin A *et al.* Bioinformatics 2013;29:15.

### Trimming / UMIs
- **fastp.** Chen S *et al.* Bioinformatics 2018;34:i884.
- **BBDuk / BBMap.** Bushnell B. LBNL 2014.
- **UMI-tools.** Smith T, Heger A, Sudbery I. Genome Res 2017;27:491.
- **Lexogen QuantSeq 3′.** Moll P *et al.* Nat Methods 2014;11 (Application Note).

### Counting + normalisation
- **SAMtools / HTSlib.** Li H *et al.* Bioinformatics 2009;25:2078.
- **HTSeq.** Anders S *et al.* Bioinformatics 2015;31:166 (HTSeq 2: Putri GH *et al.* 2022;38:2943).
- **RPKM/FPKM.** Mortazavi A *et al.* Nat Methods 2008;5:621. **TPM.** Wagner GP *et al.* Theory Biosci 2012;131:281. **Upper-quartile.** Bullard JH *et al.* BMC Bioinformatics 2010;11:94.

### QC / aggregation
- **FastQC.** Andrews S. babraham.ac.uk. **MultiQC.** Ewels P *et al.* Bioinformatics 2016;32:3047.

