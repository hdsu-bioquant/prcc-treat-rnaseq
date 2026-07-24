#!/usr/bin/env python
"""STAR ReadsPerGene.out.tab -> raw gene-count table (drops the 4 N_* summary rows).

Output columns: gene_id  unstranded  stranded_first  stranded_second
- `unstranded`      = STAR column 2 (htseq -s no)      -> the pRCC-TREAT PRIMARY column
- `stranded_first`  = STAR column 3 (htseq -s yes)     -> forward; the QuantSeq-FWD-correct column
- `stranded_second` = STAR column 4 (htseq -s reverse) -> reverse
This is the primary in-pipeline output for BOTH branches (no normalization here).
Usage: extract_star_counts.py <ReadsPerGene.out.tab> <out.tsv>
"""
import sys

inp, out = sys.argv[1], sys.argv[2]
with open(inp) as fh, open(out, "w") as o:
    o.write("gene_id\tunstranded\tstranded_first\tstranded_second\n")
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 4 or f[0].startswith("N_"):   # skip N_unmapped/N_multimapping/N_noFeature/N_ambiguous
            continue
        o.write("\t".join(f[:4]) + "\n")
