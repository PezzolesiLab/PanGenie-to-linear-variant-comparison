#!/usr/bin/env python3
import argparse
import gzip
import sys
from pathlib import Path
from collections import Counter


def open_maybe_gzip(path, mode="rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def extract_vcf_header(vcf_path):
    """
    Extract all header lines (## and #CHROM...) from an existing VCF.
    """
    header_lines = []
    with open_maybe_gzip(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                header_lines.append(line.rstrip("\n"))
            else:
                break
    if not header_lines:
        sys.exit(f"ERROR: No VCF header found in {vcf_path}")
    return header_lines


def print_minimal_header(sample_name):
    print("##fileformat=VCFv4.2")
    print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    print(
        '##INFO=<ID=ALT_PLACEHOLDER,Number=0,Type=Flag,'
        'Description="ALT allele was \'.\' and replaced with placeholder">'
    )
    print(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}")


def is_sv_marker_id(marker_id):
    """
    SV-like compact IDs:
      chrom_pos_ins_size
      chrom_pos_del_size
      chrom_pos_complex_size
    optionally with suffix:
      chrom_pos_del_size_z2
    """
    parts = marker_id.split("_")

    if len(parts) == 4 and parts[2] in {"ins", "del", "complex"} and parts[3].isdigit():
        return True

    if (
        len(parts) == 5
        and parts[2] in {"ins", "del", "complex"}
        and parts[3].isdigit()
        and parts[4].startswith("z")
        and parts[4][1:].isdigit()
    ):
        return True

    return False


def parse_non_sv_marker_id(marker_id):
    """
    Accept non-SV MarkerIDs in either form:
      chrom_pos_ref_alt
      chrom_pos_ref_alt_zN

    Returns (chrom, pos, ref, alt) or None.
    """
    parts = marker_id.split("_")

    if len(parts) == 4:
        chrom, pos, ref, alt = parts
        return chrom, pos, ref, alt

    if len(parts) == 5 and parts[4].startswith("z") and parts[4][1:].isdigit():
        chrom, pos, ref, alt, _suffix = parts
        return chrom, pos, ref, alt

    return None


def default_sidecar_paths(gwas_tsv, prefix=None):
    """
    Build default output paths for skipped/removed records.

    If prefix is provided:
      <prefix>.skipped_sv_ids.txt
      <prefix>.skipped_sv_rows.tsv
      <prefix>.removed_multisite_ids.txt
      <prefix>.removed_multisite_rows.tsv

    Otherwise use the GWAS basename.
    """
    if prefix:
        stem = prefix
    else:
        p = Path(gwas_tsv)
        name = p.name
        if name.endswith(".gz"):
            name = name[:-3]
        if name.endswith(".tsv"):
            name = name[:-4]
        elif name.endswith(".txt"):
            name = name[:-4]
        stem = str(p.with_name(name))

    return {
        "sv_ids": f"{stem}.skipped_sv_ids.txt",
        "sv_rows": f"{stem}.skipped_sv_rows.tsv",
        "multi_ids": f"{stem}.removed_multisite_ids.txt",
        "multi_rows": f"{stem}.removed_multisite_rows.tsv",
    }


def load_gwas_rows(gwas_tsv, markerid_col):
    """
    Read GWAS file into memory once so we can do a two-pass position filter.
    Returns:
      header_line, mid_idx, rows(list of raw lines), fields_list(list of split fields)
    """
    rows = []
    fields_list = []

    with open_maybe_gzip(gwas_tsv, "rt") as f:
        header = f.readline().rstrip("\n")
        if not header:
            sys.exit(f"ERROR: GWAS file {gwas_tsv} is empty.")

        cols = header.split()
        try:
            mid_idx = cols.index(markerid_col)
        except ValueError:
            sys.exit(
                f"ERROR: MarkerID column '{markerid_col}' not found.\n"
                f"Header was: {header}"
            )

        for line in f:
            if not line.strip():
                continue
            raw_line = line.rstrip("\n")
            fields = raw_line.split()
            if len(fields) <= mid_idx:
                continue
            rows.append(raw_line)
            fields_list.append(fields)

    return header, mid_idx, rows, fields_list


def load_removed_marker_ids(path):
    ids = set()
    if not path:
        return ids
    with open_maybe_gzip(path, "rt") as handle:
        first = handle.readline()
        if not first:
            return ids
        header = first.rstrip("\n").split()
        if "OriginalMarkerID" in header:
            idx = header.index("OriginalMarkerID")
            for line in handle:
                if not line.strip():
                    continue
                fields = line.rstrip("\n").split()
                if len(fields) > idx:
                    ids.add(fields[idx])
        else:
            first = first.strip()
            if first and not first.startswith("#"):
                ids.add(first.split()[0])
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                ids.add(line.split()[0])
    return ids


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Convert GWAS results with MarkerID into a single-sample presence/absence VCF.\n"
            "All genotypes are set to GT=1/1.\n\n"
            "Rules:\n"
            "  - SV-like variants are skipped entirely and written to sidecar files.\n"
            "  - Non-SV variants are kept only if their CHROM,POS is unique.\n"
            "    If a position has more than one non-SV variant, ALL variants at that position are removed\n"
            "    and written to sidecar files.\n"
            "  - Remaining records with ALT='.' use a placeholder ALT.\n"
            "Summary counts are printed to stderr."
        )
    )
    ap.add_argument("gwas_tsv", help="GWAS results file containing MarkerID column.")
    ap.add_argument("--markerid-col", default="MarkerID", help="MarkerID column name (default: MarkerID).")
    ap.add_argument("--header-vcf", default=None, help="Optional VCF to reuse header from.")
    ap.add_argument("--sample-name", default="SAMPLE", help="Sample name for output VCF.")
    ap.add_argument("--placeholder-alt", default="G", help="ALT placeholder to use (default: G).")
    ap.add_argument(
        "--sidecar-prefix",
        default=None,
        help="Optional prefix for skipped/removed sidecar files."
    )
    ap.add_argument("--skipped-sv-ids", default=None)
    ap.add_argument("--skipped-sv-rows", default=None)
    ap.add_argument("--removed-multisite-ids", default=None)
    ap.add_argument("--removed-multisite-rows", default=None)
    ap.add_argument(
        "--remove-marker-ids",
        default=None,
        help=(
            "Optional list/report of MarkerIDs to remove before VCF creation and "
            "include in removed_multisite sidecars. Tables with OriginalMarkerID "
            "are supported."
        ),
    )

    args = ap.parse_args()
    args.placeholder_alt = args.placeholder_alt.upper()

    paths = default_sidecar_paths(args.gwas_tsv, prefix=args.sidecar_prefix)
    if args.skipped_sv_ids:
        paths["sv_ids"] = args.skipped_sv_ids
    if args.skipped_sv_rows:
        paths["sv_rows"] = args.skipped_sv_rows
    if args.removed_multisite_ids:
        paths["multi_ids"] = args.removed_multisite_ids
    if args.removed_multisite_rows:
        paths["multi_rows"] = args.removed_multisite_rows

    # ---- Print VCF header ----
    if args.header_vcf:
        header_lines = extract_vcf_header(args.header_vcf)

        has_gt = any(line.startswith("##FORMAT=<ID=GT") for line in header_lines)
        has_alt_flag = any(line.startswith("##INFO=<ID=ALT_PLACEHOLDER") for line in header_lines)

        for line in header_lines:
            if line.startswith("#CHROM"):
                if not has_gt:
                    print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
                if not has_alt_flag:
                    print(
                        '##INFO=<ID=ALT_PLACEHOLDER,Number=0,Type=Flag,'
                        'Description="ALT allele was \'.\' and replaced with placeholder">'
                    )
                print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + args.sample_name)
            else:
                print(line)
    else:
        print_minimal_header(args.sample_name)

    # ---- Load GWAS once ----
    header, mid_idx, rows, fields_list = load_gwas_rows(args.gwas_tsv, args.markerid_col)
    forced_removed_ids = load_removed_marker_ids(args.remove_marker_ids)

    # ---- First pass: count non-SV variants by position ----
    pos_counts = Counter()
    n_skipped_bad_id_firstpass = 0
    n_skipped_sv_firstpass = 0

    for fields in fields_list:
        mid = fields[mid_idx]

        if mid in forced_removed_ids:
            continue

        if is_sv_marker_id(mid):
            n_skipped_sv_firstpass += 1
            continue

        parsed = parse_non_sv_marker_id(mid)
        if parsed is None:
            n_skipped_bad_id_firstpass += 1
            continue

        chrom, pos, _ref, _alt = parsed
        pos_counts[(chrom, pos)] += 1

    multi_positions = {pos for pos, count in pos_counts.items() if count > 1}

    # ---- Second pass: emit VCF + write sidecar files ----
    seen_vcf = set()
    seen_sv_ids = set()
    seen_multi_ids = set()

    n_placeholder = 0
    n_emitted = 0
    n_skipped_sv = 0
    n_skipped_bad_id = 0
    n_removed_multisite = 0
    n_forced_removed = 0

    with open(paths["sv_ids"], "w") as out_sv_ids, \
         open(paths["sv_rows"], "w") as out_sv_rows, \
         open(paths["multi_ids"], "w") as out_multi_ids, \
         open(paths["multi_rows"], "w") as out_multi_rows:

        out_sv_rows.write(header + "\n")
        out_multi_rows.write(header + "\n")

        for mid in sorted(forced_removed_ids):
            if mid not in seen_multi_ids:
                out_multi_ids.write(mid + "\n")
                seen_multi_ids.add(mid)

        for raw_line, fields in zip(rows, fields_list):
            mid = fields[mid_idx]

            # Remove IDs flagged upstream as ambiguous/colliding. These are
            # written through the same sidecars as multi-position removals so
            # downstream GWAS recovery excludes them with the existing logic.
            if mid in forced_removed_ids:
                n_removed_multisite += 1
                n_forced_removed += 1

                out_multi_rows.write(raw_line + "\n")
                continue

            # Remove SV-like records entirely, but save them
            if is_sv_marker_id(mid):
                n_skipped_sv += 1

                if mid not in seen_sv_ids:
                    out_sv_ids.write(mid + "\n")
                    seen_sv_ids.add(mid)

                out_sv_rows.write(raw_line + "\n")
                continue

            parsed = parse_non_sv_marker_id(mid)
            if parsed is None:
                n_skipped_bad_id += 1
                continue

            chrom, pos, ref, alt = parsed

            # Remove all non-SV variants at positions with >1 non-SV record
            if (chrom, pos) in multi_positions:
                n_removed_multisite += 1

                if mid not in seen_multi_ids:
                    out_multi_ids.write(mid + "\n")
                    seen_multi_ids.add(mid)

                out_multi_rows.write(raw_line + "\n")
                continue

            ref = ref.upper()
            alt = alt.upper()
            info = "."

            if (not alt) or alt == ".":
                alt = args.placeholder_alt
                info = "ALT_PLACEHOLDER"
                n_placeholder += 1

            key = (chrom, pos, ref, alt)
            if key in seen_vcf:
                continue
            seen_vcf.add(key)

            print(f"{chrom}\t{pos}\t{mid}\t{ref}\t{alt}\t.\tPASS\t{info}\tGT\t1/1")
            n_emitted += 1

    print(
        f"# Converted {n_emitted} variants; "
        f"skipped {n_skipped_sv} SV-like variants; "
        f"removed {n_removed_multisite} variants from multi-record positions; "
        f"{n_forced_removed} removed by explicit MarkerID list; "
        f"skipped {n_skipped_bad_id} malformed/nonstandard IDs; "
        f"{n_placeholder} used ALT placeholder ('{args.placeholder_alt}')",
        file=sys.stderr,
    )
    print(f"# Wrote skipped SV IDs: {paths['sv_ids']}", file=sys.stderr)
    print(f"# Wrote skipped SV rows: {paths['sv_rows']}", file=sys.stderr)
    print(f"# Wrote removed multi-position IDs: {paths['multi_ids']}", file=sys.stderr)
    print(f"# Wrote removed multi-position rows: {paths['multi_rows']}", file=sys.stderr)
    print(
        f"# First-pass counts: {len(pos_counts):,} unique non-SV positions; "
        f"{len(multi_positions):,} positions removed for having >1 non-SV variant",
        file=sys.stderr,
    )
    if forced_removed_ids:
        print(
            f"# Explicit removed MarkerIDs requested: {len(forced_removed_ids):,}",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
