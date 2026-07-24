#!/usr/bin/env python
"""Minimal GENCODE GTF -> BED12 (per transcript) for RSeQC read_distribution.py.

One BED12 line per transcript (thickStart/thickEnd from CDS where present, else
the transcript bounds). Validate against a known transcript before production use.
Usage: gtf_to_bed12.py <gencode.gtf> <out.bed12>
"""
import sys, re
from collections import defaultdict

gtf, out = sys.argv[1], sys.argv[2]
tx = {}                                   # tx_id -> dict(chrom, strand, exons[], cds[])

def attr(s, key):
    m = re.search(key + r' "([^"]+)"', s)
    return m.group(1) if m else None

with open(gtf) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        feat, chrom, start, end, strand, attrs = f[2], f[0], int(f[3]), int(f[4]), f[6], f[8]
        if feat not in ("exon", "CDS"):
            continue
        tid = attr(attrs, "transcript_id")
        if tid is None:
            continue
        t = tx.setdefault(tid, {"chrom": chrom, "strand": strand, "exon": [], "cds": []})
        t[feat.lower()].append((start - 1, end))   # BED is 0-based half-open

with open(out, "w") as o:
    for tid, t in tx.items():
        ex = sorted(t["exon"])
        if not ex:
            continue
        chrom_start = ex[0][0]
        chrom_end = ex[-1][1]
        if t["cds"]:
            cds = sorted(t["cds"])
            thick_start, thick_end = cds[0][0], cds[-1][1]
        else:
            thick_start, thick_end = chrom_end, chrom_end   # non-coding: zero-length thick
        sizes = ",".join(str(e - s) for s, e in ex) + ","
        starts = ",".join(str(s - chrom_start) for s, e in ex) + ","
        o.write("\t".join(map(str, [
            t["chrom"], chrom_start, chrom_end, tid, 0, t["strand"],
            thick_start, thick_end, 0, len(ex), sizes, starts])) + "\n")
