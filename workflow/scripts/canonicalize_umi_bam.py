#!/usr/bin/env python3
"""Canonicalize coordinate-sort ties before UMI-tools deduplication.

STAR and multithreaded ``samtools sort`` can emit the same alignments in a
slightly different order when several records share the same coordinate.
UMI-tools may need to choose one representative from otherwise equivalent
records, so a fixed random seed is only fully reproducible when those records
are presented in a deterministic order as well.

This helper expects a coordinate-sorted BAM.  It preserves the coordinate
ordering, but imposes a deterministic total order within each equal-coordinate
group.  Records with no reference coordinate (the tail of a coordinate-sorted
BAM) are name-sorted on disk and then deterministically ordered within each
query-name group, keeping memory use bounded for production-sized BAMs.

The output is an internal disposable BAM used only as the input to
``umi_tools dedup``.  The site-retained ``genomic.sorted.bam`` is not modified.
"""

from __future__ import print_function

import argparse
import os
import shutil
import tempfile

import pysam


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="coordinate-sorted input BAM")
    parser.add_argument("--output", required=True, help="canonicalized output BAM")
    parser.add_argument("--tmpdir", required=True, help="directory for temporary BAMs")
    parser.add_argument("--threads", type=int, default=1, help="threads for temporary name sort")
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("--threads must be >= 1")
    return args


def record_tiebreak_key(record):
    """Deterministic ordering for records sharing the same reference/position."""
    # Keep strand explicit before the complete SAM representation.  The SAM
    # representation supplies a total lexical order over QNAME, FLAG, CIGAR,
    # sequence/qualities and tags without inventing an analysis-specific rule.
    return (1 if record.is_reverse else 0, record.to_string())


def write_sorted_group(records, out_bam):
    if not records:
        return 0
    records.sort(key=record_tiebreak_key)
    for record in records:
        out_bam.write(record)
    return len(records)


def append_unplaced_canonically(unplaced_bam, out_bam, header, threads, workdir):
    """Append no-coordinate records in deterministic name/tie order."""
    if not os.path.exists(unplaced_bam) or os.path.getsize(unplaced_bam) == 0:
        return 0

    name_sorted = os.path.join(workdir, "unplaced.name_sorted.bam")
    sort_prefix = os.path.join(workdir, "unplaced.sorttmp")
    pysam.sort(
        "-n",
        "-@", str(threads),
        "-T", sort_prefix,
        "-o", name_sorted,
        unplaced_bam,
    )

    written = 0
    group = []
    current_name = None
    with pysam.AlignmentFile(name_sorted, "rb") as fh:
        for record in fh.fetch(until_eof=True):
            name = record.query_name or ""
            if current_name is None:
                current_name = name
            elif name != current_name:
                written += write_sorted_group(group, out_bam)
                group = []
                current_name = name
            group.append(record)
        written += write_sorted_group(group, out_bam)
    return written


def canonicalize(input_bam, output_bam, tmpdir, threads):
    os.makedirs(os.path.dirname(os.path.abspath(output_bam)), exist_ok=True)
    os.makedirs(tmpdir, exist_ok=True)

    tmp_output = output_bam + ".canonicalizing"
    if os.path.exists(tmp_output):
        os.remove(tmp_output)

    mapped_or_placed = 0
    unplaced = 0
    coordinate_groups = 0

    try:
        with tempfile.TemporaryDirectory(prefix="umi_bam_canonicalize_", dir=tmpdir) as workdir:
            unplaced_path = os.path.join(workdir, "unplaced.bam")

            with pysam.AlignmentFile(input_bam, "rb") as inp:
                header_dict = inp.header.to_dict()
                sort_order = str(header_dict.get("HD", {}).get("SO", ""))
                if sort_order and sort_order != "coordinate":
                    raise ValueError(
                        "input BAM must be coordinate sorted; header declares SO={!r}".format(sort_order)
                    )

                with pysam.AlignmentFile(tmp_output, "wb", header=inp.header) as out_bam, \
                     pysam.AlignmentFile(unplaced_path, "wb", header=inp.header) as unplaced_out:
                    current_coord = None
                    group = []
                    seen_unplaced = False

                    for record in inp.fetch(until_eof=True):
                        ref_id = int(record.reference_id)
                        ref_start = int(record.reference_start)

                        # Records without any coordinate belong at the tail of a
                        # coordinate-sorted BAM.  Sort them separately on disk so
                        # a large unmapped tail never has to be buffered in memory.
                        if ref_id < 0 or ref_start < 0:
                            seen_unplaced = True
                            unplaced_out.write(record)
                            unplaced += 1
                            continue

                        if seen_unplaced:
                            raise ValueError(
                                "input BAM is not coordinate sorted: a placed record occurs "
                                "after records without coordinates"
                            )

                        coord = (ref_id, ref_start)
                        if current_coord is None:
                            current_coord = coord
                        elif coord < current_coord:
                            raise ValueError(
                                "input BAM is not coordinate sorted: {} follows {}".format(
                                    coord, current_coord
                                )
                            )
                        elif coord != current_coord:
                            mapped_or_placed += write_sorted_group(group, out_bam)
                            coordinate_groups += 1
                            group = []
                            current_coord = coord

                        group.append(record)

                    if group:
                        mapped_or_placed += write_sorted_group(group, out_bam)
                        coordinate_groups += 1

                    # Close the unplaced writer before asking pysam/samtools to
                    # name-sort that temporary BAM.

                # Re-open the canonical output in append mode is not supported by
                # BAM writers, so rebuild it once if an unplaced tail exists.
                if unplaced:
                    placed_only = os.path.join(workdir, "placed.canonical.bam")
                    shutil.move(tmp_output, placed_only)
                    with pysam.AlignmentFile(placed_only, "rb") as placed_in, \
                         pysam.AlignmentFile(tmp_output, "wb", header=placed_in.header) as final_out:
                        for record in placed_in.fetch(until_eof=True):
                            final_out.write(record)
                        append_unplaced_canonically(
                            unplaced_path,
                            final_out,
                            placed_in.header,
                            threads,
                            workdir,
                        )

        os.replace(tmp_output, output_bam)
        # UMI-tools opens coordinate-sorted BAM inputs through pysam.fetch(),
        # which requires an index.  Build the index here so the canonical BAM
        # and its BAI are produced atomically by the same workflow step.
        index_path = output_bam + ".bai"
        if os.path.exists(index_path):
            os.remove(index_path)
        pysam.index(output_bam)
    except Exception:
        if os.path.exists(tmp_output):
            os.remove(tmp_output)
        index_path = output_bam + ".bai"
        if os.path.exists(index_path):
            os.remove(index_path)
        raise

    print(
        "Canonicalized UMI-dedup BAM input: {} placed records in {} coordinate groups; "
        "{} records without coordinates".format(
            mapped_or_placed, coordinate_groups, unplaced
        )
    )


def main():
    args = parse_args()
    canonicalize(args.input, args.output, args.tmpdir, args.threads)


if __name__ == "__main__":
    main()
