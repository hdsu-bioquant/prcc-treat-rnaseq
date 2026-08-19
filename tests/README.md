# Corley_Benchmark — testing RNA-Seq-pRCC on Corley et al. 2019

Corley et al. 2019 (*Sci Rep* 9:18895, GEO **GSE123523** / SRA **PRJNA509074**) is the
matched benchmark for the harmonization: **6 PBMC samples (3 control, 3 Poly(I:C)),
each prepared with BOTH Illumina TruSeq (full-length) AND Lexogen QuantSeq 3′** → 12
libraries. Because the *same biology* was measured with both protocols, it is the right
dataset to test whether full-length and QuantSeq outputs can be integrated.

## Matched design (the key to the sample sheet)
GEO titles `Rseq_Sx` (TruSeq) and `Qseq_Sx` (QuantSeq) share the index `Sx` = the **same
biological sample**. So `patient = S1..S6` links each TruSeq ↔ QuantSeq pair; condition is
fixed per sample (S1/S3/S5 = control, S2/S4/S6 = Poly(I:C)).

| patient | condition | full-length (TruSeq, PE 150) | QuantSeq (3′ SE 75) |
|---|---|---|---|
| S1 | control  | FL_S1 = SRR8309088 | QS_S1 = SRR8309094 |
| S2 | treated  | FL_S2 = SRR8309089 | QS_S2 = SRR8309095 |
| S3 | control  | FL_S3 = SRR8309090 | QS_S3 = SRR8309096 |
| S4 | treated  | FL_S4 = SRR8309091 | QS_S4 = SRR8309097 |
| S5 | control  | FL_S5 = SRR8309092 | QS_S5 = SRR8309098 |
| S6 | treated  | FL_S6 = SRR8309093 | QS_S6 = SRR8309099 |

Data: FASTQs at `/path/to/corley_pbmc/fastq/` (already extracted). Corley QuantSeq is **FWD, no UMI** → `has_umi: false`.

## Files here
- `samples.tsv` — the 12-sample matched sheet (absolute FASTQ paths into project 40).
- `config_corley.yaml` — pipeline config; reuses the GDC references in `RNA-Seq-pRCC/resources/gdc`; `has_umi:false`; `cross_assay.enabled:true` (produces the 3′-restricted full-length matrix, the core benchmark comparison). Results → `Corley_Benchmark/results`.

## Run (from the RNA-Seq-pRCC working dir, references reused)
```bash
conda activate snakemake-9.19.0
module load system/singularity
kinit -r 7d
cd /path/to/pRCC-RNA-Seq
snakemake -n --configfile Corley_Benchmark/config_corley.yaml          # dry-run (59 jobs)
snakemake --profile profiles/slurm --configfile Corley_Benchmark/config_corley.yaml \
          --keep-going 2>&1 | tee Corley_Benchmark/run.log
```
This is full-depth (6 PE ~92M-read + 6 SE ~30M-read STAR runs) → a real cluster job.
For a quick smoke-test first, subsample the FASTQs (e.g. first ~2M reads) into a local
`data/` and point `samples.tsv` there.

## Outputs (under `Corley_Benchmark/results/`)
- `full_length/*/*.star_gene_counts.tsv`, `quantseq/*/*.star_gene_counts.tsv` — raw STAR counts (unstranded + stranded diagnostics)
- `matrix/gene_counts_matrix.tsv` — cohort matrix (STAR unstranded, all 12)
- `matrix/three_prime_restricted_matrix.tsv` — full-length reads counted in the 3′ window (to compare against QuantSeq)
- `qc/multiqc_report.html`

Interpretation guidance is in `RNA-Seq-pRCC/docs/Pipeline_and_Corley_Explained.md`.
