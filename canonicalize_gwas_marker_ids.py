#!/usr/bin/env python3
"""Canonicalize GWAS MarkerID allele case for pipeline_v2 inputs."""

import argparse
from collections import Counter
import gzip
import sys
from pathlib import Path


def open_text(path, mode):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def canonical_marker_id(marker_id):
    parts = marker_id.split("_")
    if len(parts) < 4:
        return marker_id

    # Keep compact SV IDs unchanged; downstream type parsing expects lowercase
    # ins/del/complex tokens.
    if parts[2] in {"ins", "del", "complex"}:
        return marker_id

    # Canonical non-SV IDs: chrom_pos_REF_ALT, with optional trailing suffix.
    if len(parts) == 4 or (
        len(parts) == 5 and parts[4].startswith("z") and parts[4][1:].isdigit()
    ):
        parts[2] = parts[2].upper()
        parts[3] = parts[3].upper()
        return "_".join(parts)

    return marker_id


def header_marker_index(src):
    with open_text(src, "rt") as handle:
        header = handle.readline()
        if not header:
            return header, None, []
        cols = header.rstrip("\n").split()
        try:
            marker_idx = cols.index("MarkerID")
        except ValueError:
            raise RuntimeError(f"MarkerID column not found in {src}")
    return header, marker_idx, cols


def count_canonical_ids(src, marker_idx):
    counts = Counter()
    examples = {}
    with open_text(src, "rt") as handle:
        next(handle, None)
        for line_no, line in enumerate(handle, start=2):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split()
            if len(fields) <= marker_idx:
                continue
            old_marker = fields[marker_idx]
            new_marker = canonical_marker_id(old_marker)
            counts[new_marker] += 1
            if new_marker not in examples:
                examples[new_marker] = []
            if len(examples[new_marker]) < 5:
                examples[new_marker].append((line_no, old_marker))
    return counts, examples


def write_collision_report(report_path, collision_ids, examples):
    if not report_path:
        return
    report_path = Path(report_path)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with report_path.open("w") as handle:
        handle.write("CanonicalMarkerID\tLine\tOriginalMarkerID\n")
        for marker_id in sorted(collision_ids):
            for line_no, original in examples.get(marker_id, []):
                handle.write(f"{marker_id}\t{line_no}\t{original}\n")


def canonicalize_file(src, dst, report_path=None):
    header, marker_idx, cols = header_marker_index(src)
    if not header:
        Path(dst).write_text("")
        return 0, 0, 0

    canonical_counts, examples = count_canonical_ids(src, marker_idx)
    collision_ids = {marker_id for marker_id, count in canonical_counts.items() if count > 1}
    write_collision_report(report_path, collision_ids, examples)

    changed = 0
    removed_collisions = 0
    skipped_collisions = 0

    with open_text(src, "rt") as in_handle, open(dst, "w") as out_handle:
        in_handle.readline()

        allele_indices = [
            i
            for i, col in enumerate(cols)
            if col in {"Allele1", "Allele2", "A1", "A2", "REF", "ALT"}
        ]

        out_handle.write(header)
        for line_no, line in enumerate(in_handle, start=2):
            if not line.strip():
                out_handle.write(line)
                continue

            fields = line.rstrip("\n").split()
            if len(fields) <= marker_idx:
                out_handle.write(line)
                continue

            old_marker = fields[marker_idx]
            new_marker = canonical_marker_id(old_marker)

            if new_marker in collision_ids:
                if new_marker != old_marker:
                    skipped_collisions += 1
                removed_collisions += 1
                continue
            else:
                final_marker = new_marker

            if final_marker != old_marker:
                changed += 1
                fields[marker_idx] = final_marker

            for idx in allele_indices:
                if idx < len(fields):
                    fields[idx] = fields[idx].upper()

            out_handle.write("\t".join(fields) + "\n")

    return changed, len(collision_ids), skipped_collisions, removed_collisions


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("src", type=Path)
    parser.add_argument("dst", type=Path)
    parser.add_argument(
        "--report",
        type=Path,
        help="Optional TSV report for MarkerIDs left unchanged because uppercase canonicalization would collide.",
    )
    args = parser.parse_args()

    args.dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.dst.with_name(args.dst.name + ".tmp")

    try:
        changed, collision_groups, skipped_collisions, removed_collisions = canonicalize_file(
            args.src,
            tmp,
            args.report,
        )
    except Exception as exc:
        tmp.unlink(missing_ok=True)
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    tmp.replace(args.dst)
    msg = f"Canonicalized MarkerID allele case: {args.dst} ({changed} changed rows)"
    if collision_groups:
        msg += (
            f"; removed {removed_collisions} rows across {collision_groups} "
            "canonical MarkerID collision groups"
        )
        if skipped_collisions:
            msg += f" ({skipped_collisions} would have changed case)"
        if args.report:
            msg += f"; report: {args.report}"
    print(msg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
