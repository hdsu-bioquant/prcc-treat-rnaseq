# Synthetic smoke-test fixture

This directory contains a deliberately tiny, fully synthetic RNA-seq fixture for the
pRCC-RNA-Seq pipeline. It is intended to live in Git and to run quickly during development
and at partner sites. It tests **pipeline routing and the production output contract**, not
biological realism or GDC-reference compatibility; those belong in `tests/real/`.

## What the fixture covers

The miniature reference (`reference/synthetic.fa` + `synthetic.gtf`) contains three
protein-coding genes and supports spliced STAR alignment. Four libraries exercise every
current core route:

| library_id | layout | assay | UMI | purpose |
|---|---|---|---|---|
| `FL_noUMI` | paired | full length | no | baseline FL path + splice alignment |
| `FL_UMI` | paired | full length | 6 nt at R1 start | proves UMI logic is assay-independent |
| `QS_noUMI` | single | QuantSeq | no | QuantSeq trimming/alignment path |
| `QS_UMI` | single | QuantSeq | 6 nt at R1 start | UMI extraction + dedup + QuantSeq trimming |

The four libraries deliberately map to only two `sample_id` values, testing that multiple
libraries may belong to one biological sample. The extra `batch` column in `samples.tsv`
demonstrates that additional metadata columns are accepted without controlling workflow
routing.

The UMI fixtures deliberately contain all three important deduplication cases:

1. same mapping position + same UMI -> collapse;
2. same mapping position + different UMI -> retain;
3. same UMI + different mapping position -> retain.

QuantSeq reads include terminal poly(A); each QS library also contains one deliberately
unmappable read and one read that becomes too short after trimming.

## Expected values

`expected/` contains the independent golden numerical expectations:

- `gene_lengths.tsv`: union-exon lengths and gene metadata;
- `raw_gene_counts.tsv`: exact STAR **unstranded** raw counts;
- `umi_dedup_gene_counts.tsv`: exact UMI molecule counts;
- `expected_summary.tsv`: input/post-trim/assigned/molecule totals;
- `read_manifest.tsv`: read-level blueprint for debugging.

The fixture FASTQs and reference are tiny enough to commit directly to Git.
`checksums.sha256` verifies the **input fixture** before every run.

## Running the test

From anywhere, run:

```bash
bash /path/to/pRCC-RNA-Seq/tests/synthetic/run_test.sh
```

The script resolves the repository root automatically and then:

1. records basic execution information in a timestamped local log;
2. verifies the committed synthetic fixture using `checksums.sha256`;
3. removes any previous `tests/synthetic/output/`;
4. removes and rebuilds `tests/synthetic/reference/star_index/` and generated `gene_lengths.tsv`;
5. runs the normal Snakemake workflow with `tests/synthetic/config.yaml`;
6. validates the resulting **production-style output structure** and exact expected values.

A successful run ends with:

```text
Synthetic smoke test PASSED.
```

## What a successful synthetic run looks like

The test uses the same output contract as a normal analysis:

```text
tests/synthetic/output/
├── results/                         # portable canonical result package
│   ├── libraries/<library_id>/
│   │   ├── gene_expression.tsv
│   │   └── qc_metrics.tsv
│   ├── matrices/
│   │   ├── raw_gene_counts.tsv
│   │   └── umi_molecule_counts.tsv
│   ├── qc/
│   │   ├── qc_metrics.tsv
│   │   └── multiqc_report.html
│   └── run/
│       ├── libraries.tsv
│       ├── config.yaml
│       ├── provenance.yaml
│       ├── software_versions.tsv
│       ├── references.tsv
│       ├── manifest.tsv
│       ├── checksums.sha256
│       └── validation_checksums.sha256
├── restricted/                    # site-retained BAMs/detailed QC/full local metadata
└── intermediate/                  # disposable processing artefacts
```

`validate_results.py` checks the exact raw counts and UMI molecule counts, the canonical
`gene_expression.tsv` schema (including FL normalization vs QS `NA` values), selected stable
QC metrics, the library-level MultiQC summary, complete manifest-driven software-version
provenance, restricted BAM/FastQC presence, run metadata, and both generated checksum files.
The default mixed-level MultiQC General Statistics table is intentionally suppressed; FastQC
R1/R2 diagnostics remain in their dedicated report sections.

The synthetic DAG adds only one provenance job for software versions. It records the software
subset declared by `workflow/config/software_versions.yaml`, writes both
`results/run/software_versions.tsv` and MultiQC's `*_mqc_versions.yml`, and records the actual
Snakemake controller version. It does not launch per-container version-probe jobs.

MultiQC 1.21 may still print `custom_content | prcc_rnaseq_qc: Found 1 samples (table)` for
the single custom-content table file. That message counts the parsed custom-content object, not
the library rows; the rendered Library QC Summary and validator require all four library rows.

`results/run/checksums.sha256` is for **package/transfer integrity**. The separate
`validation_checksums.sha256` contains only deliberately deterministic canonical products
and is intended for cross-installation comparison. MultiQC HTML and timestamped provenance
are excluded from validation hashes.

Once the output contract is frozen for a consortium release, a reference
`validation_checksums.sha256` can be committed under `expected/`; partners can then compare
their synthetic run byte-for-byte for this deterministic subset in addition to the current
numerical validation.

## Partner-site qualification

For an initial installation check, a partner can send the maintainers:

- `tests/synthetic/output/results/` (the same portable package used for real runs); and
- the timestamped `tests/synthetic/logs/run_test_*.log` if requested.

This provides evidence that all core assay/layout/UMI paths and the data-delivery contract
work at that site. A separate real-data harmonization run under `tests/real/` is intended to
validate the production GDC reference/resources and normal site-specific execution profile.

## Local execution logs

Every invocation mirrors stdout/stderr to:

```text
tests/synthetic/logs/run_test_YYYYMMDD_HHMMSS.log
```

These logs are Git-ignored. They may contain host names, usernames, filesystem paths, and
other infrastructure details emitted by Snakemake/tools, so review them before public
sharing.

## Git tracking policy

Tracked fixture/source files include the synthetic FASTQs, FASTA/GTF, sample sheet,
configuration, expected tables, generator, validator, and `run_test.sh`.

Generated locally and **not** tracked:

- `output/`;
- `reference/star_index/`;
- `reference/gene_lengths.tsv`;
- `logs/`.

## Manual debugging

Equivalent manual execution from the repository root is:

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --configfile tests/synthetic/config.yaml \
  --cores 2 \
  --software-deployment-method apptainer \
  --printshellcmds

python tests/synthetic/validate_results.py
```

For a true clean run, remove `tests/synthetic/output/`, `tests/synthetic/reference/star_index/`, and
`tests/synthetic/reference/gene_lengths.tsv` first; `run_test.sh` does this automatically.

## Regeneration

```bash
python tests/synthetic/generate_synthetic_data.py
```

Gzipped FASTQs are deterministic (`mtime=0`). Regeneration also rewrites the synthetic
config/sample sheet/golden fixture tables and refreshes `checksums.sha256`.
