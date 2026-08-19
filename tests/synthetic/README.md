# Synthetic smoke-test fixture

This directory contains a deliberately tiny, fully synthetic RNA-seq fixture for the
pRCC-TREAT pipeline. It is intended to live in Git and to run quickly at every development
iteration/CI run. It is **not** intended to test biological realism or GDC-reference
compatibility; those belong in `tests/real/`.

## Why a miniature reference is included

The synthetic reads are generated against `reference/synthetic.fa` and
`reference/synthetic.gtf`. This makes exact expected gene counts known in advance and
avoids requiring the ~25-GB GDC STAR index for every smoke test. The real-data suite will
continue to test the genuine GRCh38.d1.vd1/GENCODE v36 path.

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
metrics because metadata/order/timestamps can legitimately differ.

## Running the smoke test

From the repository root, run the workflow with the synthetic configuration using the
same Snakemake/container runtime used for the production pipeline. For example, with
Apptainer enabled:

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile tests/synthetic/config.yaml \
  --cores 2 --software-deployment-method apptainer
```

Then validate the canonical numerical outputs:

```bash
python tests/synthetic/validate_results.py
```

The validator requires exact agreement for the raw gene-count and UMI molecule-count
matrices and also checks that MultiQC was produced.

## Regeneration

Run `python tests/synthetic/generate_synthetic_data.py` from the repository root. Gzipped
FASTQs are written deterministically (`mtime=0`). `checksums.sha256` records the committed
fixture bytes.
