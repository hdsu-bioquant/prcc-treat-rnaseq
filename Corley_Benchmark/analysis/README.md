# Corley benchmark — reproducible analysis

This folder makes the Corley interpretation **re-runnable**. Until now the numbers in
`RNA-Seq-pRCC/docs/Pipeline_and_Corley_Explained.md` §5.2 were computed ad-hoc; this
script regenerates them from the pipeline outputs with explicit, documented choices.

## Run

```bash
conda run -n snakemake-9.19.0 python corley_analysis.py     # needs pandas + numpy (no scipy)
```

Inputs (all pipeline outputs): `results/matrix/gene_counts_matrix.tsv`,
`results/matrix/three_prime_restricted_matrix.tsv`,
`results/qc/multiqc_report_data/multiqc_star.txt`, `samples.tsv`, and
`RNA-Seq-pRCC/resources/gdc/gencode.v36.gene_lengths.tsv` (gene_id ↔ gene_name).

## Design
6 patients (S1–S6), each sequenced with **both** protocols. `control` = S1/S3/S5,
`treated` (Poly(I:C)) = S2/S4/S6. Full-length ↔ QuantSeq matched by `patient`. Batch constant.

## Outputs (`tables/`)

| File | What it shows |
|---|---|
| `qc_mapping_summary.tsv` | STAR uniquely-mapped % per library (final/pass-2) |
| `spearman_correlation_12x12.tsv` | full 12×12 Spearman matrix (rank + Pearson) |
| `correlation_summary.tsv` | within- vs cross-protocol means + nearest-neighbour test |
| `isg_induction.tsv` | per-ISG control→treated fold change, **both** protocols |
| `isg_signature_per_sample.tsv` | mean ISG log2CPM per library |
| `three_prime_restriction.tsv` | FL↔QS agreement: full-gene vs 3′-restricted |

## Key results (reproducible)

- **Mapping:** full-length ≈ 94.4–94.9 %, QuantSeq ≈ 88.2–91.0 % uniquely mapped.
- **Protocol is the dominant axis:** within-FL Spearman **0.944**, within-QS **0.910**,
  cross-protocol matched **0.827** — and **every** library's nearest correlate is the
  *same* protocol (**12/12**). Integration therefore needs protocol correction.
- **Both protocols recover the Poly(I:C) response, concordantly:** 30/30 canonical ISGs
  detected; **median induction FL 18.2× vs QuantSeq 17.9×**; ISG signature 5.99 → 10.24 log2CPM
  (control → treated, full-length). Top ISGs agree across protocols (IFI27 134×/141×,
  IFIT1 66×/81×, RSAD2 63×/72×, CXCL10 53×/65×).
- **3′-restriction hurts:** FL↔QS agreement falls from **0.827** (full-gene) to **0.69**
  (3′-restricted), Δ −0.14 → the 3′-restricted secondary is *not* the route to integration.

## Method notes
- **Spearman** computed as Pearson-on-ranks (avoids the scipy gap in the snakemake env).
- **CPM** = raw / library-size × 1e6; correlations on expressed genes (row sum > 0, 43 779 genes).
- **ISG panel** (30 canonical type-I interferon-stimulated genes) mapped by symbol → GENCODE gene_id.
- Uses the pipeline's **unstranded** primary counts (the archival matrix) for both protocols.

## Note on this copy
This is the copy under `40_pRCC-TREAT/pipelines/pRCC-RNA-Seq/`. Paths in
`corley_analysis.py` are derived from the script's own location, so it works
wherever the pipeline lives. Differential expression is **not** included in this
copy (it was run in the 38_Pipelines working copy).
