# Synthetic smoke-test fixture

This directory contains a deliberately tiny, fully synthetic RNA-seq fixture for the
pRCC-TREAT pipeline. It is intended to live in Git and to run quickly at every development
iteration and at partner sites. It is **not** intended to test biological realism or
GDC-reference compatibility; those belong in `tests/real/`.

## Why a miniature reference is included

The synthetic reads are generated against `reference/synthetic.fa` and
`reference/synthetic.gtf`. This makes exact expected gene counts known in advance and
avoids requiring the production GDC STAR index for every smoke test. The real-data suite
will continue to test the genuine GRCh38/GENCODE path.

The reference contains three non-overlapping protein-coding test genes:

- `SYN_GENE_A` / `SYN_A`: plus strand, two exons
- `SYN_GENE_B` / `SYN_B`: minus strand, two exons
- `SYN_GENE_C` / `SYN_C`: plus strand, one exon

Some full-length fragments cross A/B splice junctions, so the test exercises spliced STAR
alignment rather than only contiguous genomic mapping.

## Four route-covering libraries

| sample | layout | assay | UMI | purpose |
|---|---|---|---|---|
| `FL_noUMI` | PE | full length | no | baseline FL path + splice alignment |
| `FL_UMI` | PE | full length | 6 nt at R1 start | proves UMI logic is not coupled to QS |
| `QS_noUMI` | SE | 3' | no | QuantSeq trimming/alignment path |
| `QS_UMI` | SE | 3' | 6 nt at R1 start | Lexogen-like UMI extraction + dedup |

The UMI fixtures deliberately contain all three cases required to test deduplication:

1. same mapping position + same UMI -> PCR duplicates should collapse;
2. same mapping position + different UMI -> distinct molecules should remain;
3. same UMI + different mapping position -> distinct molecules should remain.

QuantSeq reads additionally carry a terminal poly(A) segment. Each QS sample contains one
unmappable read and one deliberately over-trimmed read that should fail `minlength=20`.

## Exact expectations

- `expected/gene_lengths.tsv`: expected union-exon gene lengths from the miniature GTF.
- `expected/raw_gene_counts.tsv`: expected STAR **unstranded raw** gene counts.
- `expected/umi_dedup_gene_counts.tsv`: expected molecule-level counts after UMI dedup.
- `expected/expected_summary.tsv`: expected read/fragment totals, trimming loss, and unmapped controls.
- `expected/read_manifest.tsv`: read-level blueprint for diagnosing failed assertions.

Do not use whole-file checksums for every pipeline output. Canonical count tables can be
compared exactly; BAM/log/HTML outputs should be checked structurally or by selected parsed
metrics because metadata, ordering, compression, and timestamps can legitimately differ.

## Running the smoke test

From the repository root, run:

```bash
bash tests/synthetic/run_test.sh
```

The script resolves the repository root from its own location, so it can also be invoked
from another working directory as long as you provide a valid path to `run_test.sh`.

`run_test.sh` performs a clean end-to-end test:

1. records basic provenance (date, host, Git commit/state, Snakemake/Python/container runtime versions);
2. verifies the committed fixture against `checksums.sha256`;
3. removes any previous `tests/synthetic/results/` directory;
4. removes and rebuilds `tests/synthetic/reference/star_index/`;
5. runs Snakemake with `tests/synthetic/config.yaml` using Apptainer;
6. runs `validate_results.py` and requires exact agreement with the golden raw-count and
   UMI-deduplicated molecule-count matrices, plus a MultiQC report.

A successful run ends with:

```text
Synthetic smoke test PASSED.
```

### Local run logs

Every invocation writes the complete terminal output (stdout and stderr) to a timestamped
file while still displaying it interactively:

```text
tests/synthetic/logs/run_test_YYYYMMDD_HHMMSS.log
```

The `logs/` directory is local and Git-ignored. These files are intended for troubleshooting
and can be shared with the pipeline maintainers when a partner installation fails.

**Note:** execution logs can contain host names, usernames, filesystem paths, project names,
and other local infrastructure information emitted by Snakemake or external tools. Review a
log before posting it publicly or attaching it to a public GitHub issue.

## Manual debugging

If necessary, the two main steps can also be run manually from the repository root:

```bash
snakemake   --snakefile workflow/Snakefile   --configfile tests/synthetic/config.yaml   --cores 2   --software-deployment-method apptainer   --printshellcmds

python tests/synthetic/validate_results.py
```

For a true clean smoke test, remove both `tests/synthetic/results/` and
`tests/synthetic/reference/star_index/` first; `run_test.sh` does this automatically.

## Regeneration

Run:

```bash
python tests/synthetic/generate_synthetic_data.py
```

from the repository root. Gzipped FASTQs are written deterministically (`mtime=0`).
`checksums.sha256` records the committed fixture bytes and excludes generated `results/`,
`reference/star_index/`, and `logs/` content.
