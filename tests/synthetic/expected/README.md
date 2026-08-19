# Expected outputs

`raw_gene_counts.tsv` and `umi_dedup_gene_counts.tsv` are the main golden numerical
outputs. They should be compared by parsed table content rather than by relying on the
formatting of upstream STAR/HTSeq files.

The molecule counts assume exact-coordinate UMI-aware deduplication with UMIs extracted
from the first six bases of R1. Deliberately distinct UMIs at the same coordinate have
large Hamming distances, so adjacency-based UMI-tools methods should not merge them.
