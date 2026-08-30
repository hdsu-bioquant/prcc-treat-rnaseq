# Consortium IO definition — sample sheet

## Contract

One row represents one sequencing library. `library_id` is the workflow/output identifier. `sample_id` is the biological sample identifier and may repeat when multiple libraries represent the same biological sample.

Always-required columns:

```text
library_id  sample_id  assay  layout  strandedness  fq1  fq2  has_umi
```

If at least one library has `has_umi=true`, the sheet must also contain:

```text
umi_pattern  umi_location  umi_discard_bases
```

## Supported consortium combinations

- `assay=full_length`, `layout=paired`
- `assay=quantseq`, `layout=single`

Supported strandedness values are `unstranded`, `forward`, and `reverse` and must reflect the library-preparation protocol.

## Identifiers

`library_id` must be unique. `library_id` and `sample_id` follow:

```text
^[A-Za-z0-9][A-Za-z0-9._-]*$
```

## FASTQs

`fq1` is required. `fq2` is required for paired libraries and absent/`-` for single-end libraries. Supported suffixes are `.fastq`, `.fastq.gz`, `.fq`, and `.fq.gz`.

The same FASTQ must not be assigned to more than one library, including symlink-equivalent paths.

## UMI fields

UMI handling is library-specific and assay-independent.

- `has_umi`: `true` or `false`.
- `umi_pattern`: one or more `N` characters describing a fixed contiguous 5′ UMI.
- `umi_location`: `read1_start` or `read2_start`; R2 requires paired-end layout.
- `umi_discard_bases`: non-negative integer number of immediately adjacent bases removed after UMI extraction but not used as part of the molecular identifier.

## Additional columns

Additional metadata columns may be present and are retained in the loaded run table, but they do not currently alter workflow routing. The original sheet is retained under `restricted/run/`; the portable library manifest contains only standardized technical fields.

## Validation

Validate before execution:

```bash
python workflow/scripts/sample_sheet.py /path/to/run/samples.tsv
```
