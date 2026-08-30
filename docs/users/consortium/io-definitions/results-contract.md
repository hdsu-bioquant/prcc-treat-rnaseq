# Consortium IO definition — results contract

| Field | Value |
|---|---|
| Document ID | IO-03 |
| Status | Draft |
| Document version | 0.1 |
| Owner | pRCC-TREAT bioinformatics maintainers |
| Applicable pipeline release | `development` |
| Last revised | 2026-08-30 |

## Output contract version

The maintained output contract is version 1.

```text
output/
  results/
  restricted/
  intermediate/
```

Only `results/` is the standard portable consortium package.

## Canonical per-library expression table

```text
results/libraries/<library_id>/gene_expression.tsv
```

Columns, in maintained order:

```text
gene_id gene_name gene_type gene_length unstranded stranded_first stranded_second fpkm fpkm_uq tpm umi_molecule_count
```

Semantics:

- STAR gene-count columns are preserved as unstranded/first/second strand counts;
- full-length libraries receive GDC-compatible FPKM, FPKM-UQ, and TPM values;
- QuantSeq normalized columns are `NA`;
- UMI libraries receive molecule counts;
- non-UMI libraries have `NA` molecule counts.

## Run-level matrices

```text
results/matrices/raw_gene_counts.tsv
```

is produced for supported runs.

When UMI libraries are present:

```text
results/matrices/umi_molecule_counts.tsv
```

is also produced.

## Boundaries

`restricted/` is not part of the portable contract and remains local unless separately approved. `intermediate/` is workflow working state and is not delivered.

A pipeline release can change without changing output-contract version. The output-contract version changes only when the maintained portable interface changes incompatibly.
