# Synthetic smoke-test fixture

This directory contains a deliberately tiny, fully synthetic RNA-seq fixture for the
pRCC-TREAT-RNA-Seq pipeline. It is intended to live in Git and to run quickly during development
and at partner sites. It tests **pipeline routing and the production output contract**, not
biological realism or GDC-reference compatibility; those belong in `tests/real/`.

## What the fixture covers

The miniature reference (`reference/synthetic.fa` + `synthetic.gtf`) contains three
protein-coding genes and supports spliced STAR alignment. Four libraries exercise every
current core route:

| library_id | layout | assay | UMI | purpose |
|---|---|---|---|---|
| `FL_noUMI` | paired | full length | no | baseline FL path + splice alignment |
| `FL_UMI` | paired | full length | 8 nt at R2 start, discard 0 | non-Lexogen UMI architecture + dedup |
| `QS_noUMI` | single | QuantSeq | no | QuantSeq trimming/alignment path |
| `QS_UMI` | single | QuantSeq | 6 nt at R1 start + 4 nt discard | Lexogen-like UMI + spacer extraction, dedup + QuantSeq trimming |

The four libraries deliberately map to only two `sample_id` values, testing that multiple
libraries may belong to one biological sample. The extra `batch` column in `samples.tsv`
demonstrates that additional metadata columns are accepted without controlling workflow
routing.

The UMI fixtures deliberately use **different extraction specifications**. `FL_UMI` carries an
8-base UMI at the start of R2 with no adjacent discard bases. `QS_UMI` carries a 6-base UMI
at the start of R1 followed by the literal four-base synthetic spacer `TATA`; the sample sheet
describes that spacer only as `umi_discard_bases=4`, so the workflow must not depend on its
sequence identity. This mirrors the real Lexogen structure that motivated the UMI refactor while
proving that the implementation is not hard-coded to Lexogen.

The UMI fixtures also deliberately contain all three important deduplication cases:

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
- `read_manifest.tsv`: read-level blueprint for debugging;
- `validation_checksums.sha256`: frozen maintainer-approved hashes for the 14 deterministic
  canonical run products used for cross-installation qualification.

The fixture FASTQs and reference are tiny enough to commit directly to Git.
`checksums.sha256` verifies the **input/golden fixture** before every run. It deliberately does
**not** include `expected/validation_checksums.sha256`, because that file is the separately
maintained deterministic workflow-output baseline. Updating an output baseline therefore does not
require regenerating the synthetic fixture or its integrity manifest.

## Running the test

From anywhere, run:

```bash
bash /path/to/prcc-treat-rnaseq/tests/synthetic/run_test.sh
```

The script resolves the repository root automatically and then:

1. records basic execution information in a timestamped local log;
2. verifies the committed synthetic fixture using `checksums.sha256`;
3. removes any previous `tests/synthetic/output/`;
4. removes and rebuilds `tests/synthetic/reference/star_index/` and generated `gene_lengths.tsv`;
5. runs the normal Snakemake workflow with `tests/synthetic/config.yaml`, retaining `temp()`
   intermediates with `--notemp` so UMI extraction can be inspected directly;
6. validates the exact raw-to-extracted UMI FASTQ transformation before checking alignment/counts;
7. validates the resulting **production-style output structure** and exact expected values;
8. verifies the generated deterministic output hashes against the frozen reference baseline in
   `expected/validation_checksums.sha256`.

A successful run ends with:

```text
Synthetic smoke test PASSED.
```


### Container/software upgrades

The synthetic validator also checks the maintained UMI-tools declaration. The current qualified
version is **UMI-tools 1.1.6**, with deduplication using the workflow-owned `--random-seed=1` after
deterministic canonicalization of equal-coordinate BAM ties.
After changing a pinned container, refresh the local image with `bash containers/pull_images.sh` before
running this test. First run the normal smoke test against the maintained deterministic baseline. If
the generated validation manifest is unchanged, no additional reproducibility run is required.

### Maintainer baseline-update mode

When an intentional change is expected to alter deterministic canonical files, first run the normal
smoke test and inspect any baseline mismatch. If that workflow completed successfully and failed
only on comparison with the old frozen baseline, validate the existing output again with:

```bash
python tests/synthetic/validate_results.py --skip-frozen-baseline
```

This rechecks fixture semantics, exact UMI extraction, biological expectations, the output contract,
package checksums, and the newly generated deterministic validation manifest while skipping only
comparison with `expected/validation_checksums.sha256`. A second workflow execution is unnecessary
for a reviewed metadata/provenance-only mismatch.

Preview or accept an intentional candidate baseline with:

```bash
bash tests/maintainers/update_validation_baseline.sh synthetic
bash tests/maintainers/update_validation_baseline.sh synthetic --apply
```

Do not use the helper to update `tests/synthetic/checksums.sha256`; that separate manifest protects
the version-controlled input/golden fixture and intentionally excludes the deterministic output
baseline. Metadata/provenance-only baseline differences may be accepted after review. If
deterministic computational products changed, two fresh clean runs of the finalized implementation
must produce byte-identical generated validation manifests before `--apply` is used.
See `docs/maintainers/qualification-baselines.md` for the complete maintainer procedure. Normal
qualification runs should omit `--skip-frozen-baseline`.

### Execution-environment independence

The synthetic qualification deliberately does **not** use the production workstation or
SLURM profile templates. `run_test.sh` disables inherited Snakemake global/workflow profiles
and runs the fixed tiny DAG locally with two cores. This keeps the qualification procedure
identical across partner sites and prevents site-specific scheduler/resource settings from
changing what is being tested.

On an HPC system where computation on the login node is not permitted, obtain a compute-node
allocation (interactive or batch, according to local policy) and run the same
`bash tests/synthetic/run_test.sh` command there. The site-specific production profile is
validated separately by `tests/real/`, where execution is intentionally production-like.

The rule-level `mem_mb` values remain the production workflow's scheduling declarations. The
synthetic test does not set a global `mem_mb` scheduling limit, so those declarations do not
reserve that amount of physical RAM during direct local execution. Actual memory use is that
of the miniature synthetic reference and tiny inputs. A 32-GB workstation is therefore a
reasonable environment for this qualification, but passing it does **not** qualify that
machine for full human GDC-scale processing.

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

`validate_results.py` first verifies the two declared UMI architectures against the committed raw
FASTQs and then checks the retained UMI-extracted FASTQs byte-semantically: `FL_UMI` must remove
exactly the first 8 bases from R2 while leaving R1 sequence/quality unchanged, and `QS_UMI` must
remove exactly 6 UMI + 4 spacer bases from R1. It also verifies that the extracted UMI is propagated
into processed read names for downstream UMI-aware deduplication.

It then checks the exact raw counts and UMI molecule counts, the canonical
`gene_expression.tsv` schema (including FL normalization vs QS `NA` values), selected stable
QC metrics, the production UMI-extraction conformance QC for both synthetic UMI architectures,
the library-level MultiQC summary, complete manifest-driven software-version provenance,
restricted BAM/FastQC presence, run metadata, both generated checksum files, and the frozen
cross-installation validation baseline.
The default mixed-level MultiQC General Statistics table is intentionally suppressed; FastQC
R1/R2 diagnostics remain in their dedicated report sections.

The synthetic DAG adds only one provenance job for software versions. It records the software
subset declared by `workflow/config/software_versions.yaml`, writes both
`results/run/software_versions.tsv` and MultiQC's `*_mqc_versions.yml`, and records the actual
Snakemake controller version. It does not launch per-container version-probe jobs.

MultiQC 1.21 may still print `custom_content | prcc_treat_rnaseq_qc: Found 1 samples (table)` for
the single custom-content table file. That message counts the parsed custom-content object, not
the library rows; the rendered Library QC Summary and validator require all four library rows.

`results/run/checksums.sha256` is for **package/transfer integrity**. The separate
`validation_checksums.sha256` contains only deliberately deterministic canonical products
and is intended for cross-installation comparison. MultiQC HTML and timestamped provenance
are excluded from validation hashes.

The maintained `expected/validation_checksums.sha256` makes a successful synthetic run demonstrate
both numerical/schema correctness and byte-for-byte agreement of the 14 deterministic canonical
products with the expected baseline. It must not be silently regenerated: intentional
metadata/provenance differences may be accepted after review, while a new computational output
baseline requires the reproducibility procedure described above.

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
  --profile none \
  --workflow-profile none \
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
config/sample sheet/golden fixture tables and refreshes `checksums.sha256`. That fixture manifest
does **not** include or rewrite `expected/validation_checksums.sha256`; the deterministic
workflow-output baseline is maintained separately with the maintainer baseline-update procedure.
