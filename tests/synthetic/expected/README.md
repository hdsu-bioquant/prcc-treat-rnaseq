# Expected outputs

`raw_gene_counts.tsv` and `umi_dedup_gene_counts.tsv` are the main golden numerical
outputs. They should be compared by parsed table content rather than by relying on the
formatting of upstream STAR/HTSeq files.

The molecule counts assume exact-coordinate UMI-aware deduplication after the fixture-specific
UMI extraction step. `FL_UMI` uses an 8-nt UMI at the start of R2 with no additional discard
bases. `QS_UMI` uses a 6-nt UMI at the start of R1 followed by four bases that must be discarded.
Deliberately distinct UMIs at the same coordinate have large Hamming distances, so adjacency-based
UMI-tools methods should not merge them.

`validation_checksums.sha256` is the maintainer-approved byte-level reference for the
14 deterministic canonical files emitted by a successful synthetic run. It is intentionally
**not generated or refreshed by `generate_synthetic_data.py`**. Updating this file is a
maintenance action that should only happen after intentional deterministic output changes have
been reviewed and reproduced in clean runs.

