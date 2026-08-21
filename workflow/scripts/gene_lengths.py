#!/usr/bin/env python3
"""Build the compact gene-information table used for expression normalization.

Gene length is the union length of all annotated exons for a gene.  The table
also carries gene type/name and chromosome because the current GDC FPKM-UQ
calculation defines its upper-quartile denominator on autosomal protein-coding
genes.

Usage: gene_lengths.py <annotation.gtf> <out.tsv>
Output columns: gene_id  gene_length  gene_type  gene_name  chromosome
"""

import re
import sys
from collections import defaultdict

gtf, out = sys.argv[1], sys.argv[2]
exons = defaultdict(list)  # gene_id -> [(start, end), ...]
meta = {}                  # gene_id -> (gene_type, gene_name, chromosome)


def attr(text, key):
    m = re.search(key + r' "([^"]+)"', text)
    return m.group(1) if m else "NA"


with open(gtf) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue
        chromosome = fields[0]
        feature = fields[2]
        start = int(fields[3])
        end = int(fields[4])
        attrs = fields[8]
        gid = attr(attrs, "gene_id")

        if feature == "gene":
            meta[gid] = (
                attr(attrs, "gene_type"),
                attr(attrs, "gene_name"),
                chromosome,
            )
        elif feature == "exon":
            exons[gid].append((start, end))
            # Some valid GTFs omit explicit gene features.  Preserve enough
            # metadata from exons to keep the normalization table usable.
            meta.setdefault(
                gid,
                (attr(attrs, "gene_type"), attr(attrs, "gene_name"), chromosome),
            )


def union_len(intervals):
    intervals = sorted(intervals)
    total = 0
    current_start = None
    current_end = None
    for start, end in intervals:
        if current_start is None:
            current_start, current_end = start, end
        elif start <= current_end + 1:
            current_end = max(current_end, end)
        else:
            total += current_end - current_start + 1
            current_start, current_end = start, end
    if current_start is not None:
        total += current_end - current_start + 1
    return total


with open(out, "w") as fh:
    fh.write("gene_id\tgene_length\tgene_type\tgene_name\tchromosome\n")
    for gid in sorted(exons):
        gene_type, gene_name, chromosome = meta.get(gid, ("NA", "NA", "NA"))
        fh.write(
            f"{gid}\t{union_len(exons[gid])}\t{gene_type}\t{gene_name}\t{chromosome}\n"
        )
