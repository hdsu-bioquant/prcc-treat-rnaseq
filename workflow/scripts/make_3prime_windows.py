#!/usr/bin/env python
"""Emit a GTF of the 3'-most W bp of EXONIC sequence per gene (for HTSeq-count).

Used by the Branch-A 3'-restricted secondary: counting full-length reads only in this
3' window emulates QuantSeq 3'-tag counting, removing the transcript-length bias so
full-length and QuantSeq gene counts become comparable. Strand-aware.
Emits BOTH a GTF (HTSeq counting features: `exon` features carrying gene_id;
HTSeq-count groups by gene_id via -t exon -i gene_id) AND a BED of the same windows
(used by samtools view -L to pre-filter reads before HTSeq).
NOTE: validate on a few known genes before production.
Usage: make_3prime_windows.py <gencode.gtf> <out.gtf> <out.bed> <window_bp>
"""
import sys, re
from collections import defaultdict

gtf, out, bed, W = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
exons = defaultdict(list)
meta = {}

def attr(s, k):
    m = re.search(k + r' "([^"]+)"', s)
    return m.group(1) if m else None

with open(gtf) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        gid = attr(f[8], "gene_id")
        if gid is None:
            continue
        exons[gid].append((int(f[3]), int(f[4])))
        meta[gid] = (f[0], f[6])          # chrom, strand

def merge(iv):
    iv = sorted(iv); m = []
    for s, e in iv:
        if m and s <= m[-1][1] + 1:
            m[-1] = (m[-1][0], max(m[-1][1], e))
        else:
            m.append((s, e))
    return m

with open(out, "w") as o, open(bed, "w") as b:
    for gid, iv in exons.items():
        chrom, strand = meta[gid]
        m = merge(iv)
        need = W
        rows = []
        order = m if strand == "-" else m[::-1]     # walk from the 3' end
        for s, e in order:
            if need <= 0:
                break
            L = e - s + 1
            if L <= need:
                rows.append((s, e)); need -= L
            elif strand == "-":                      # 3' = low coordinate -> take from s upward
                rows.append((s, s + need - 1)); need = 0
            else:                                    # 3' = high coordinate -> take from e downward
                rows.append((e - need + 1, e)); need = 0
        for s, e in rows:
            attrs = 'gene_id "%s"; transcript_id "%s";' % (gid, gid)
            o.write("\t".join([chrom, "prcc3p", "exon", str(s), str(e), ".", strand, ".", attrs]) + "\n")
            b.write("\t".join([chrom, str(s - 1), str(e), gid, ".", strand]) + "\n")   # BED: 0-based half-open
