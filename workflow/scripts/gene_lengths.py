#!/usr/bin/env python
"""Union-exon (non-overlapping) gene length + biotype from a GENCODE GTF.

Used as the "Gene Length (L)" term in the GDC FPKM / FPKM-UQ / TPM formulas and to
identify protein-coding genes for the FPKM/FPKM-UQ denominators.
Ref: https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/
Usage: gene_lengths.py <gencode.gtf> <out.tsv>
Output columns: gene_id  gene_length  gene_type  gene_name
"""
import sys, re
from collections import defaultdict

gtf, out = sys.argv[1], sys.argv[2]
exons = defaultdict(list)            # gene_id -> [(start, end), ...]
meta  = {}                           # gene_id -> (gene_type, gene_name)

def attr(s, key):
    m = re.search(key + r' "([^"]+)"', s)
    return m.group(1) if m else "NA"

with open(gtf) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        feat, start, end, attrs = f[2], int(f[3]), int(f[4]), f[8]
        gid = attr(attrs, "gene_id")
        if feat == "gene":
            meta[gid] = (attr(attrs, "gene_type"), attr(attrs, "gene_name"))
        elif feat == "exon":
            exons[gid].append((start, end))

def union_len(intervals):
    intervals.sort()
    total, cs, ce = 0, None, None
    for s, e in intervals:
        if cs is None:
            cs, ce = s, e
        elif s <= ce + 1:
            ce = max(ce, e)
        else:
            total += ce - cs + 1
            cs, ce = s, e
    if cs is not None:
        total += ce - cs + 1
    return total

with open(out, "w") as o:
    o.write("gene_id\tgene_length\tgene_type\tgene_name\n")
    for gid, iv in exons.items():
        gt, gn = meta.get(gid, ("NA", "NA"))
        o.write("%s\t%d\t%s\t%s\n" % (gid, union_len(iv), gt, gn))
