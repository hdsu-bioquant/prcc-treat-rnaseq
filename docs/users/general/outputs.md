# General outputs

Each run uses the configured output root:

```text
output/
  results/
  restricted/
  intermediate/
```

`results/` is the portable result package. `restricted/` can contain original inputs/metadata and optional outputs that should not be assumed portable. `intermediate/` contains workflow working files.

## Canonical expression output

For each library:

```text
results/libraries/<library_id>/gene_expression.tsv
```

Columns:

```text
gene_id gene_name gene_type gene_length unstranded stranded_first stranded_second fpkm fpkm_uq tpm umi_molecule_count
```

Full-length libraries receive GDC-compatible FPKM, FPKM-UQ, and TPM values. QuantSeq normalized columns are `NA`. UMI libraries receive molecule counts; non-UMI libraries have `NA` molecule counts.

Run-level matrices include `results/matrices/raw_gene_counts.tsv` and, when UMI libraries are present, `results/matrices/umi_molecule_counts.tsv`.

QC and provenance are stored beneath `results/qc/` and `results/run/`. The exact consortium delivery contract is defined separately under `docs/users/consortium/io-definitions/`.
