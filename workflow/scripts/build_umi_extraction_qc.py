#!/usr/bin/env python3
"""Verify that UMI extraction conforms to the declared pipeline-native specification.

The check is deliberately generic and assay-independent.  It compares a deterministic
prefix of the raw FASTQ(s) with the corresponding UMI-tools output FASTQ(s), verifying:

* the UMI-bearing read lost exactly ``umi_length + umi_discard_bases`` 5' bases;
* the non-UMI mate of a paired library was sequence/quality preserving;
* retained read names carry the expected extracted UMI using UMI-tools' read-id form.

The script reports a small machine-readable QC row.  Exact transform/tag mismatches are
pipeline-integrity errors and cause a non-zero exit.  Sampled retention is reported but
is not thresholded because short/malformed raw reads can legitimately fail extraction.

Only the first ``--max-records`` raw records are sampled (default: 10,000), keeping the
check cheap even for production FASTQs.
"""

from __future__ import annotations

import argparse
import gzip
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

import pandas as pd


@dataclass(frozen=True)
class FastqRecord:
    header: str
    sequence: str
    quality: str


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--raw-r1", required=True)
    p.add_argument("--raw-r2", default="-")
    p.add_argument("--extracted-r1", required=True)
    p.add_argument("--extracted-r2", default="-")
    p.add_argument("--umi-pattern", required=True)
    p.add_argument("--umi-location", choices=("read1_start", "read2_start"), required=True)
    p.add_argument("--umi-discard-bases", type=int, required=True)
    p.add_argument("--max-records", type=int, default=10_000)
    p.add_argument("--output", required=True)
    args = p.parse_args()
    if not re.fullmatch(r"N+", args.umi_pattern.upper()):
        p.error("--umi-pattern currently supports only fixed-length N+ specifications")
    if args.umi_discard_bases < 0:
        p.error("--umi-discard-bases must be non-negative")
    if args.max_records <= 0:
        p.error("--max-records must be positive")
    if args.umi_location == "read2_start" and args.raw_r2 == "-":
        p.error("read2_start requires paired raw FASTQs")
    if (args.raw_r2 == "-") != (args.extracted_r2 == "-"):
        p.error("raw/extracted R2 must either both be supplied or both be '-'")
    return args


def open_fastq(path: str) -> TextIO:
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def read_record(fh: TextIO, label: str) -> FastqRecord | None:
    header = fh.readline()
    if not header:
        return None
    sequence = fh.readline()
    plus = fh.readline()
    quality = fh.readline()
    if not sequence or not plus or not quality:
        raise ValueError(f"truncated FASTQ record in {label}")
    header = header.rstrip("\r\n")
    sequence = sequence.rstrip("\r\n")
    plus = plus.rstrip("\r\n")
    quality = quality.rstrip("\r\n")
    if not header.startswith("@") or not plus.startswith("+"):
        raise ValueError(f"malformed FASTQ record in {label}: {header!r}")
    if len(sequence) != len(quality):
        raise ValueError(f"sequence/quality length mismatch in {label}: {header!r}")
    return FastqRecord(header=header[1:], sequence=sequence, quality=quality)


def token(header: str) -> str:
    return header.split()[0]


def normalize_pair_suffix(read_id: str) -> str:
    return re.sub(r"/[12]$", "", read_id)


def raw_id(header: str) -> str:
    return normalize_pair_suffix(token(header))


def processed_identity(header: str, umi_length: int) -> tuple[str, str | None]:
    """Return (base read id, appended UMI or None) from an extracted read header."""
    t = token(header)

    # Be tolerant of a trailing /1 or /2 if an upstream naming convention keeps it
    # after the UMI.  UMI-tools normally appends ``_<UMI>`` as the last read-id word.
    trailing_pair = ""
    m = re.search(r"(/[12])$", t)
    if m:
        trailing_pair = m.group(1)
        t = t[: -len(trailing_pair)]

    suffix_len = umi_length + 1
    if len(t) >= suffix_len and t[-suffix_len] == "_":
        observed_umi = t[-umi_length:]
        base = t[:-suffix_len]
        return normalize_pair_suffix(base), observed_umi

    # If the expected UMI suffix is absent, preserve the best available base ID so
    # the sampled record can still be counted as retained but tag-nonconforming.
    return normalize_pair_suffix(t + trailing_pair), None


def load_raw_sample(args: argparse.Namespace) -> dict[str, dict[str, object]]:
    paired = args.raw_r2 != "-"
    umi_length = len(args.umi_pattern)
    remove_n = umi_length + args.umi_discard_bases
    sample: dict[str, dict[str, object]] = {}

    with open_fastq(args.raw_r1) as r1_fh:
        r2_ctx = open_fastq(args.raw_r2) if paired else None
        try:
            for _ in range(args.max_records):
                r1 = read_record(r1_fh, args.raw_r1)
                r2 = read_record(r2_ctx, args.raw_r2) if paired else None
                if r1 is None:
                    if paired and r2 is not None:
                        raise ValueError("raw paired FASTQs have different record counts")
                    break
                if paired and r2 is None:
                    raise ValueError("raw paired FASTQs have different record counts")

                rid = raw_id(r1.header)
                if paired and raw_id(r2.header) != rid:
                    raise ValueError(
                        f"raw paired FASTQ IDs differ: {r1.header!r} vs {r2.header!r}"
                    )
                if rid in sample:
                    raise ValueError(f"duplicate raw read ID in QC sample: {rid}")

                primary = r1 if args.umi_location == "read1_start" else r2
                if len(primary.sequence) < remove_n:
                    expected_umi = primary.sequence[:umi_length]
                    eligible = False
                else:
                    expected_umi = primary.sequence[:umi_length]
                    eligible = True

                sample[rid] = {
                    "r1": r1,
                    "r2": r2,
                    "expected_umi": expected_umi,
                    "eligible": eligible,
                    "remove_n": remove_n,
                }
        finally:
            if r2_ctx is not None:
                r2_ctx.close()

    if not sample:
        raise ValueError("raw FASTQ contains no records")
    return sample


def compare_output(args: argparse.Namespace, sample: dict[str, dict[str, object]]) -> tuple[int, int, int]:
    paired = args.extracted_r2 != "-"
    umi_length = len(args.umi_pattern)
    found: set[str] = set()
    transform_matches = 0
    tag_matches = 0

    # Every retained record originating in the first N raw records must occur within
    # the first N extracted records because extraction preserves order while only
    # potentially dropping records.
    with open_fastq(args.extracted_r1) as e1_fh:
        e2_ctx = open_fastq(args.extracted_r2) if paired else None
        try:
            for _ in range(len(sample)):
                e1 = read_record(e1_fh, args.extracted_r1)
                e2 = read_record(e2_ctx, args.extracted_r2) if paired else None
                if e1 is None:
                    if paired and e2 is not None:
                        raise ValueError("extracted paired FASTQs have different record counts")
                    break
                if paired and e2 is None:
                    raise ValueError("extracted paired FASTQs have different record counts")

                id1, tag1 = processed_identity(e1.header, umi_length)
                if paired:
                    id2, tag2 = processed_identity(e2.header, umi_length)
                    if id1 != id2:
                        raise ValueError(
                            f"extracted paired FASTQ IDs differ: {e1.header!r} vs {e2.header!r}"
                        )
                else:
                    id2, tag2 = id1, tag1

                rid = id1
                if rid not in sample:
                    # This is normally the first record beyond the sampled raw prefix
                    # when one or more sampled reads were dropped. Keep scanning until
                    # the N-record output bound is reached.
                    continue
                if rid in found:
                    raise ValueError(f"duplicate extracted read ID in QC sample: {rid}")
                found.add(rid)

                raw = sample[rid]
                expected_umi = str(raw["expected_umi"])
                r1: FastqRecord = raw["r1"]  # type: ignore[assignment]
                r2: FastqRecord | None = raw["r2"]  # type: ignore[assignment]
                remove_n = int(raw["remove_n"])

                tags_ok = tag1 == expected_umi and (not paired or tag2 == expected_umi)
                if tags_ok:
                    tag_matches += 1

                if args.umi_location == "read1_start":
                    r1_ok = (
                        e1.sequence == r1.sequence[remove_n:]
                        and e1.quality == r1.quality[remove_n:]
                    )
                    r2_ok = True if not paired else (
                        e2.sequence == r2.sequence and e2.quality == r2.quality
                    )
                else:
                    r1_ok = e1.sequence == r1.sequence and e1.quality == r1.quality
                    r2_ok = (
                        e2.sequence == r2.sequence[remove_n:]
                        and e2.quality == r2.quality[remove_n:]
                    )

                if r1_ok and r2_ok:
                    transform_matches += 1
        finally:
            if e2_ctx is not None:
                e2_ctx.close()

    return len(found), transform_matches, tag_matches


def pct(n: int, d: int):
    return (n / d * 100.0) if d else pd.NA


def main() -> None:
    args = parse_args()
    sample = load_raw_sample(args)
    retained, transform_matches, tag_matches = compare_output(args, sample)

    row = {
        "umi_length": len(args.umi_pattern),
        "umi_location": args.umi_location,
        "umi_discard_bases": args.umi_discard_bases,
        "umi_extract_qc_records": len(sample),
        "umi_extract_qc_retained_percent": pct(retained, len(sample)),
        "umi_extract_transform_match_percent": pct(transform_matches, retained),
        "umi_extract_tag_match_percent": pct(tag_matches, retained),
    }

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([row]).to_csv(out, sep="\t", index=False, na_rep="NA", float_format="%.6g")

    if retained == 0:
        raise SystemExit("UMI extraction QC found no retained records in the deterministic sample")
    if transform_matches != retained:
        raise SystemExit(
            "UMI extraction QC failed: retained FASTQ sequence/quality transformation does not "
            f"match the declared UMI specification for {retained - transform_matches}/{retained} records"
        )
    if tag_matches != retained:
        raise SystemExit(
            "UMI extraction QC failed: extracted read-name UMI tag does not match the raw prefix for "
            f"{retained - tag_matches}/{retained} records"
        )


if __name__ == "__main__":
    main()
