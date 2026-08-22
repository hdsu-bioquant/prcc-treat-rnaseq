# Expected outputs

`raw_gene_counts.tsv` and `umi_dedup_gene_counts.tsv` are the main golden numerical
outputs. They should be compared by parsed table content rather than by relying on the
formatting of upstream STAR/HTSeq files.

The molecule counts assume exact-coordinate UMI-aware deduplication with UMIs extracted
from the first six bases of R1. Deliberately distinct UMIs at the same coordinate have
large Hamming distances, so adjacency-based UMI-tools methods should not merge them.

`validation_checksums.sha256` is the maintainer-approved byte-level reference for the
14 deterministic canonical files emitted by a successful synthetic run. It is intentionally
**not generated or refreshed by `generate_synthetic_data.py`**. Updating this file is a
release/maintenance action that should only happen after the deterministic output contract
changes have been reviewed and reproduced in clean runs.

