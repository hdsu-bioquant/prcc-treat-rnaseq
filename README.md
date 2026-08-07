# RNA-Seq-pRCC pipeline — User guide

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

The sample sheet (`config/samples.tsv`) is **tab-separated WITH a header row**. One row per **library**:

```
sample   assay   fq1   fq2   patient   batch   strandedness   [condition]
```


## 3. Pipeline configuration

Everything is driven by `config/config.yaml`.

### 3.1 Main settings

```yaml
samples: config/samples.tsv            
results: results                       
tmpdir:  results/tmp                  

reference:
  dir: resources/gdc                   
  download_references: true            
  genome_fasta: GRCh38.d1.vd1.fa
  gtf:          gencode.v36.annotation.gtf
  star_index:   star-2.7.5c_GRCh38.d1.vd1_gencode.v36
  sjdb_overhang: 100                   
```

### 3.2 Assay branch settings

```yaml
full_length:
  trim_adapters:    false              
  count_column:     unstranded         
  compute_fpkm_tpm: true               

quantseq:
  has_umi:      true                   
  umi_pattern:  "NNNNNN"               
  strandedness: forward                
  bbduk_polyA:  true                  
```



### 3.3 Tunable parameters

```yaml
star:
  threads: 8
  gdc_params: >-                       # the VERBATIM GDC DR32+ STAR recipe — do NOT edit in order to preserve GDC conformity
    --twopassMode Basic --quantMode TranscriptomeSAM GeneCounts
    --outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNoverLmax 0.1
    --alignSJoverhangMin 8 --chimSegmentMin 15  [... full GDC block ...]
```

---

## 4. Pipeline execution

Ideal start status: inside `tmux`, Apptainer/Singularity available, `snakemake-9.19.0` active, working
directory at the pipeline root, and references (§0.5) + images (§0.6) present.

### 4.1 Dry-run

```bash
snakemake -n --configfile config/config.yaml
```

### 4.2 Launching the pipeline

Depending on the **infrastructure** of your working environment, the **snakemake workflow** can be edited to use a **profile** personalised to your architecture: **HPC cluster** or **local working station** (see the steps below):

**On an HPC cluster (SLURM):**

```bash
snakemake --profile profiles/slurm \
          --configfile config/config.yaml \
          --rerun-incomplete --keep-going
```


> **⚙ Adapt the SLURM profile to your cluster** — edit `profiles/slurm/config.yaml`


**On a local workstation:**

```bash
snakemake --profile profiles/local \
          --configfile config/config.yaml \
          --rerun-incomplete --keep-going
```

---

## 5. Outputs & Interpretation

Under `results/`:

```
qc/multiqc_report.html                                aggregated QC (FastQC + STAR)
qc/fastqc/<sample>.done                               per-sample FastQC marker (+ html/zip)

full_length/<s>/
  ├── <s>.Aligned.sortedByCoord.bam(.bai)             coord-sorted genome BAM
  ├── <s>.ReadsPerGene.out.tab                        raw STAR GeneCounts (3 strand columns)
  ├── <s>.star_gene_counts.tsv                        ★ PRIMARY per-sample counts (cleaned)
  └── <s>.augmented_star_gene_counts.tsv              FPKM / FPKM-UQ / TPM   (if compute_fpkm_tpm)

quantseq/<s>/
  ├── <s>.star_gene_counts.tsv                        ★ PRIMARY QuantSeq counts (non-dedup, uniform basis)
  ├── <s>.dedup.bam                                   UMI-deduplicated BAM (secondary path)
  └── <s>.dedup_htseq_counts.tsv                      SECONDARY UMI-dedup HTSeq counts

matrix/
  ├── gene_counts_matrix.tsv                          ★★ PRIMARY cohort matrix (STAR unstranded, BOTH branches)
  └── quantseq_umidedup_matrix.tsv                    SECONDARY QuantSeq UMI-dedup matrix (if has_umi)
```

## 6. Optional modules

### 6.1 FPKM / FPKM-UQ / TPM (GDC augmented output)

```yaml
full_length:
  compute_fpkm_tpm: true    
```
Output:
`full_length/<s>/<s>.augmented_star_gene_counts.tsv`.

### 6.2 UMI-dedup secondary (QuantSeq)

```yaml
quantseq:
  has_umi:      true
  umi_pattern:  "NNNNNN"     # 6 bp UMI at the 5' of Read 1
  strandedness: forward      # HTSeq -s: forward→yes | reverse→reverse | unstranded→no
```
Outputs:
`quantseq/<s>/<s>.dedup_htseq_counts.tsv`, `matrix/quantseq_umidedup_matrix.tsv`.

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
| `quay.io/biocontainers/umi_tools:1.1.4--py39hf95cd2a_2` | UMI extract + dedup (QuantSeq) |
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

## References

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

