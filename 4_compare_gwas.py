#!/usr/bin/env python3
import argparse
import gzip
import sys
import math
from pathlib import Path
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# -------------------------
# Helpers
# -------------------------
def open_maybe_gzip(path, mode="rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def canonical_chrom(chrom: str) -> str:
    return chrom[3:] if chrom.lower().startswith("chr") else chrom


def parse_marker_id(marker_id):
    """
    MarkerID formats expected:
      Non-SV-like: chrom_pos_ref_alt
                   chrom_pos_ref_alt_zN
      SV-like:     chrom_pos_ins_size
                   chrom_pos_del_size
                   chrom_pos_complex_size
                   optionally with suffix _zN

    Returns (chrom, pos)
    """
    parts = marker_id.split("_")
    if len(parts) < 2:
        raise ValueError
    chrom = canonical_chrom(parts[0])
    pos = int(parts[1])
    return chrom, pos


def classify_variant(marker_id: str) -> str:
    """
    Classify MarkerID into SNV, INDEL, SV, or UNKNOWN.
    """
    parts = marker_id.split("_")

    if len(parts) < 4:
        return "UNKNOWN"

    # SV-like compact IDs
    if parts[2] in {"ins", "del", "complex"}:
        if len(parts) == 4 and parts[3].isdigit():
            return "SV"
        if (
            len(parts) == 5
            and parts[3].isdigit()
            and parts[4].startswith("z")
            and parts[4][1:].isdigit()
        ):
            return "SV"
        return "UNKNOWN"

    # Non-SV IDs
    if len(parts) == 4:
        ref = parts[2]
        alt = parts[3]
    elif len(parts) == 5 and parts[4].startswith("z") and parts[4][1:].isdigit():
        ref = parts[2]
        alt = parts[3]
    else:
        return "UNKNOWN"

    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    return "INDEL"


def safe_neg_log10(p):
    if p <= 0:
        p = 1e-300
    elif p > 1:
        p = 1.0
    return -math.log10(p)


def read_id_list(path):
    ids = set()
    with open_maybe_gzip(path, "rt") as f:
        for line in f:
            line = line.strip()
            if line:
                ids.add(line)
    return ids


# -------------------------
# Load GWAS for unique recovery
# -------------------------
def load_gwas_file(path):
    """
    Load GWAS file into:
      gwas_by_chrom[chrom][pos] = list of GWAS lines
      gwas_by_id[MarkerID] = original GWAS line
    """
    path = Path(path)
    with open_maybe_gzip(path, "rt") as f:
        header = f.readline().rstrip("\n")
        if not header:
            sys.exit(f"ERROR: File {path} is empty.")

        cols = header.split()
        try:
            mid_idx = cols.index("MarkerID")
        except ValueError:
            sys.exit(f"ERROR: MarkerID column not found in {path}")

        gwas_by_chrom = defaultdict(lambda: defaultdict(list))
        gwas_by_id = {}

        for line in f:
            if not line.strip():
                continue
            raw = line.rstrip("\n")
            fields = raw.split()
            if len(fields) <= mid_idx:
                continue
            mid = fields[mid_idx]
            try:
                chrom, pos = parse_marker_id(mid)
            except ValueError:
                continue
            gwas_by_chrom[chrom][pos].append(raw)
            gwas_by_id[mid] = raw

    return header, mid_idx, gwas_by_chrom, gwas_by_id


# -------------------------
# Load GWAS for shared metrics creation
# -------------------------
def load_gwas_effects_by_pos(path):
    """
    Mirrors the shared beta/pvalue script logic.
    Returns:
      header, gwas[(chrom,pos)] = list of recs
    rec fields:
      MarkerID, BETA, SE, p, AF
    """
    path = Path(path)
    with open_maybe_gzip(path, "rt") as f:
        header = f.readline().rstrip("\n")
        if not header:
            sys.exit(f"ERROR: {path} empty")

        cols = header.split()
        req = ["MarkerID", "BETA", "SE", "p.value", "AF_Allele2"]
        try:
            idx = {c: cols.index(c) for c in req}
        except ValueError as e:
            sys.exit(f"ERROR: missing required column in {path}: {e}")

        gwas = defaultdict(list)
        for line in f:
            if not line.strip():
                continue
            fields = line.split()
            try:
                mid = fields[idx["MarkerID"]]
                chrom, pos = parse_marker_id(mid)
                rec = {
                    "MarkerID": mid,
                    "BETA": float(fields[idx["BETA"]]),
                    "SE": float(fields[idx["SE"]]),
                    "p": float(fields[idx["p.value"]]),
                    "AF": float(fields[idx["AF_Allele2"]]),
                }
            except Exception:
                continue
            gwas[(chrom, pos)].append(rec)

    return header, gwas


# -------------------------
# Read hap.py lists as positions
# -------------------------
def read_happy_positions(path):
    """
    hap.py IDs -> set of (chrom, pos)
    """
    pos_set = set()
    with open_maybe_gzip(path, "rt") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                chrom, pos = parse_marker_id(line)
                pos_set.add((chrom, pos))
            except ValueError:
                continue
    return pos_set


# -------------------------
# Window joins
# -------------------------
def collect_gwas_rows_with_window_filtered(gwas_by_chrom, hap_pos_set, window, marker_idx, excluded_ids=None):
    excluded_ids = excluded_ids or set()
    collected = set()

    for chrom, pos in hap_pos_set:
        chrom_dict = gwas_by_chrom.get(chrom)
        if not chrom_dict:
            continue
        for p in range(pos - window, pos + window + 1):
            if p in chrom_dict:
                for line in chrom_dict[p]:
                    fields = line.split()
                    if len(fields) <= marker_idx:
                        continue
                    mid = fields[marker_idx]
                    if mid in excluded_ids:
                        continue
                    collected.add(line)

    return list(collected)


def collect_with_window(gwas_by_pos, hap_pos, window):
    """
    Mirrors the shared beta/pvalue script logic.
    Dedupes by exact MarkerID.
    """
    collected = {}
    for chrom, pos in hap_pos:
        for p in range(pos - window, pos + window + 1):
            for rec in gwas_by_pos.get((chrom, p), []):
                collected[rec["MarkerID"]] = rec
    return list(collected.values())


def add_back_ids(gwas_by_id, ids_to_add, excluded_ids=None):
    excluded_ids = excluded_ids or set()
    added = []
    missing = 0

    for mid in ids_to_add:
        if mid in excluded_ids:
            continue
        line = gwas_by_id.get(mid)
        if line is not None:
            added.append(line)
        else:
            missing += 1

    return added, missing


# -------------------------
# Map helpers
# -------------------------
def lines_to_id_map(lines, marker_idx):
    out = {}
    for line in lines:
        fields = line.split()
        if len(fields) <= marker_idx:
            continue
        out[fields[marker_idx]] = line
    return out


def shared_records_to_id_map(records):
    return {r["MarkerID"]: r for r in records}


# -------------------------
# Shared metrics writers
# -------------------------
def write_shared_preconcordance_outputs(shared_map, prefix):
    """
    shared_map[mid] = (r_lin, r_pan)
    Writes:
      prefix.shared_gwas_ids.pre_concordance.txt
      prefix.shared_variants_metrics.pre_concordance.tsv
    """
    ids_path = f"{prefix}.shared_gwas_ids.pre_concordance.txt"
    with open(ids_path, "w") as f:
        for mid in sorted(shared_map):
            f.write(mid + "\n")

    out_tsv = f"{prefix}.shared_variants_metrics.pre_concordance.tsv"
    with open(out_tsv, "w") as f:
        f.write(
            "MarkerID\t"
            "BETA_pan\tSE_pan\tp_pan\tAF_pan\t"
            "BETA_lin\tSE_lin\tp_lin\tAF_lin\t"
            "delta_beta\tdelta_logp\tdelta_af\n"
        )
        for mid in sorted(shared_map):
            r_lin, r_pan = shared_map[mid]
            db = r_pan["BETA"] - r_lin["BETA"]
            dlp = safe_neg_log10(r_pan["p"]) - safe_neg_log10(r_lin["p"])
            daf = r_pan["AF"] - r_lin["AF"]
            f.write(
                f"{mid}\t"
                f"{r_pan['BETA']}\t{r_pan['SE']}\t{r_pan['p']}\t{r_pan['AF']}\t"
                f"{r_lin['BETA']}\t{r_lin['SE']}\t{r_lin['p']}\t{r_lin['AF']}\t"
                f"{db}\t{dlp}\t{daf}\n"
            )

    print(f"Wrote: {ids_path}")
    print(f"Wrote: {out_tsv}")


# -------------------------
# AF + variant type extraction
# -------------------------
def extract_af_and_type(header, marker_idx, gwas_lines, af_col="AF_Allele2"):
    cols = header.split()
    try:
        af_idx = cols.index(af_col)
    except ValueError:
        sys.exit(f"ERROR: {af_col} not found in GWAS header.")

    results = []
    unknown_lines = []

    for line in gwas_lines:
        fields = line.split()

        try:
            af = float(fields[af_idx])
        except (ValueError, IndexError):
            continue

        if not (0.0 <= af <= 1.0):
            continue

        try:
            mid = fields[marker_idx]
        except IndexError:
            unknown_lines.append(line)
            continue

        vtype = classify_variant(mid)

        if vtype == "UNKNOWN":
            unknown_lines.append(line)
            continue

        results.append((af, vtype))

    return results, unknown_lines


# -------------------------
# AF binning
# -------------------------
def bin_af_values_by_type(af_type_list):
    step = 0.05
    edges = [i * step for i in range(0, 20)] + [1.000001]
    labels = [f"{int(edges[i]*100)}–{int(edges[i+1]*100)}%" for i in range(len(edges)-1)]

    counts = {
        "SNV": [0] * len(labels),
        "INDEL": [0] * len(labels),
        "SV": [0] * len(labels),
    }

    for af, vt in af_type_list:
        for i in range(len(labels)):
            if edges[i] <= af < edges[i + 1]:
                counts[vt][i] += 1
                break

    return labels, counts


def write_bin_type_summary(labels, counts, out_path, label):
    totals = [sum(vals) for vals in zip(*counts.values())]
    grand = sum(totals)

    with open(out_path, "w") as f:
        f.write("File\tBin\tVariantType\tCount\tBinTotal\tPercentOfBin\tPercentOfTotal\n")
        for i, binlab in enumerate(labels):
            for vt in ("SNV", "INDEL", "SV"):
                c = counts[vt][i]
                bt = totals[i]
                f.write(
                    f"{label}\t{binlab}\t{vt}\t{c}\t{bt}\t"
                    f"{(100 * c / bt if bt else 0):.4f}\t"
                    f"{(100 * c / grand if grand else 0):.4f}\n"
                )
        for vt in ("SNV", "INDEL", "SV"):
            c = sum(counts[vt])
            f.write(
                f"{label}\tAll\t{vt}\t{c}\t{grand}\t"
                f"{(100 * c / grand if grand else 0):.4f}\t"
                f"{(100 * c / grand if grand else 0):.4f}\n"
            )


# -------------------------
# Plotting for final unique sets
# -------------------------
def plot_af_histograms(
    header1, idx1, gwas_lines1,
    header2, idx2, gwas_lines2,
    name1, name2, prefix
):
    af1, unknown1 = extract_af_and_type(header1, idx1, gwas_lines1)
    af2, unknown2 = extract_af_and_type(header2, idx2, gwas_lines2)

    with open(f"{prefix}.unknown_file1_rows.tsv", "w") as f:
        f.write(header1 + "\n")
        for line in unknown1:
            f.write(line + "\n")

    with open(f"{prefix}.unknown_file2_rows.tsv", "w") as f:
        f.write(header2 + "\n")
        for line in unknown2:
            f.write(line + "\n")

    print(f"Unknown rows in {name1}: {len(unknown1):,}")
    print(f"Unknown rows in {name2}: {len(unknown2):,}")

    pdf = f"{prefix}.AF_unique_histograms_by_type.pdf"
    with PdfPages(pdf) as out:
        for af, name, tag in [(af1, name1, "file1"), (af2, name2, "file2")]:
            if not af:
                continue

            labels, counts = bin_af_values_by_type(af)
            write_bin_type_summary(
                labels,
                counts,
                f"{prefix}.unique_{tag}_AF_bin_type_summary.tsv",
                name
            )

            fig, ax = plt.subplots(figsize=(10, 4))
            bottom = [0] * len(labels)

            for vt in ("SNV", "INDEL", "SV"):
                ax.bar(range(len(labels)), counts[vt], bottom=bottom, label=vt)
                bottom = [b + c for b, c in zip(bottom, counts[vt])]

            ax.set_xticks(range(len(labels)))
            ax.set_xticklabels(labels, rotation=90, fontsize=6)
            ax.set_ylabel("Variant count")
            ax.set_title(f"AF_Allele2 (unique to {name})")
            ax.legend()
            fig.tight_layout()
            out.savefig(fig)
            plt.close(fig)

    print(f"Wrote: {pdf}")


# -------------------------
# Main
# -------------------------
def main():
    ap = argparse.ArgumentParser(
        description=(
            "Recover Pangenie-unique, Linear-unique, and shared pre-concordance sets; "
            "then reconcile exact MarkerID overlaps so the three final sets are disjoint."
        )
    )
    ap.add_argument("file1", help="Pangenie GWAS results")
    ap.add_argument("file2", help="Linear GWAS results")

    ap.add_argument("--pangenie-unique", required=True, help="hap.py FP list")
    ap.add_argument("--linear-unique", required=True, help="hap.py FN list")
    ap.add_argument("--shared-pre", required=True, help="hap.py shared pre-concordance list")

    ap.add_argument("--window", type=int, default=10, help="Window size in bp (default ±10)")
    ap.add_argument("--prefix", default=None)

    ap.add_argument(
        "--pangenie-skipped-sv-ids",
        default="1_create_happy_pangenie.skipped_sv_ids.txt",
        help="Skipped SV IDs to reinclude for pangenie unique recovery"
    )
    ap.add_argument(
        "--pangenie-removed-multisite-ids",
        default="1_create_happy_pangenie.removed_multisite_ids.txt",
        help="Removed multisite IDs to exclude from pangenie recovery"
    )
    ap.add_argument(
        "--linear-removed-multisite-ids",
        default="1_create_happy_linear.removed_multisite_ids.txt",
        help="Removed multisite IDs to exclude from linear recovery"
    )

    args = ap.parse_args()

    # Load GWAS in two forms:
    #   1) line-based recovery
    #   2) effect-based shared metrics creation
    header1, idx1, gwas1, gwas1_by_id = load_gwas_file(args.file1)
    header2, idx2, gwas2, gwas2_by_id = load_gwas_file(args.file2)

    _, pan_effects = load_gwas_effects_by_pos(args.file1)
    _, lin_effects = load_gwas_effects_by_pos(args.file2)

    pan_pos = read_happy_positions(args.pangenie_unique)
    lin_pos = read_happy_positions(args.linear_unique)
    shared_pos = read_happy_positions(args.shared_pre)

    pan_skipped_sv_ids = read_id_list(args.pangenie_skipped_sv_ids)
    pan_removed_multisite_ids = read_id_list(args.pangenie_removed_multisite_ids)
    lin_removed_multisite_ids = read_id_list(args.linear_removed_multisite_ids)

    base1 = Path(args.file1).name
    base2 = Path(args.file2).name
    prefix = args.prefix or f"{base1}__vs__{base2}"

    # --------------------------------
    # A) Recover pangenie-unique exactly as before
    # --------------------------------
    pan_lines_window = collect_gwas_rows_with_window_filtered(
        gwas1, pan_pos, args.window, idx1, excluded_ids=pan_removed_multisite_ids
    )

    pan_sv_lines, pan_sv_missing = add_back_ids(
        gwas1_by_id,
        pan_skipped_sv_ids,
        excluded_ids=pan_removed_multisite_ids
    )

    pan_lines_initial = list(lines_to_id_map(pan_lines_window + pan_sv_lines, idx1).values())
    pan_map = lines_to_id_map(pan_lines_initial, idx1)

    # --------------------------------
    # B) Recover linear-unique exactly as before
    # --------------------------------
    lin_lines_initial = collect_gwas_rows_with_window_filtered(
        gwas2, lin_pos, args.window, idx2, excluded_ids=lin_removed_multisite_ids
    )
    lin_map = lines_to_id_map(lin_lines_initial, idx2)

    # --------------------------------
    # C) Build pre-concordance shared exactly like the shared script
    # --------------------------------
    lin_rows_shared = collect_with_window(lin_effects, shared_pos, args.window)

    pan_by_id = {}
    for recs in pan_effects.values():
        for r in recs:
            pan_by_id[r["MarkerID"]] = r

    shared_initial = {}
    for r_lin in lin_rows_shared:
        mid = r_lin["MarkerID"]
        if mid in pan_by_id:
            shared_initial[mid] = (r_lin, pan_by_id[mid])

    # --------------------------------
    # Before-filter stats
    # --------------------------------
    pan_ids_before = set(pan_map)
    lin_ids_before = set(lin_map)
    shared_ids_before = set(shared_initial)

    print("Before overlap filtering:")
    print(f"  Pangenie-unique recovered IDs: {len(pan_ids_before):,}")
    print(f"  Linear-unique recovered IDs:   {len(lin_ids_before):,}")
    print(f"  Shared pre-concordance IDs:    {len(shared_ids_before):,}")

    # overlap diagnostics before filtering
    pan_shared_overlap = pan_ids_before & shared_ids_before
    lin_shared_overlap = lin_ids_before & shared_ids_before
    pan_lin_overlap = pan_ids_before & lin_ids_before
    triple_overlap = pan_ids_before & lin_ids_before & shared_ids_before

    print("Overlap counts before filtering:")
    print(f"  pan ∩ shared:   {len(pan_shared_overlap):,}")
    print(f"  lin ∩ shared:   {len(lin_shared_overlap):,}")
    print(f"  pan ∩ lin:      {len(pan_lin_overlap):,}")
    print(f"  pan ∩ lin ∩ shared: {len(triple_overlap):,}")

    # --------------------------------
    # D) Reconcile overlaps by exact MarkerID
    # --------------------------------
    shared_final = dict(shared_initial)

    # Rule 1 + 2:
    # if in shared and a unique set, remove from unique and keep in shared
    for mid in pan_shared_overlap:
        if mid in pan_map:
            del pan_map[mid]

    for mid in lin_shared_overlap:
        if mid in lin_map:
            del lin_map[mid]

    # Rule 3:
    # if in both unique sets, remove from both and move to shared
    # if not already in shared, create shared record from both GWAS effect tables
    pan_lin_only_overlap = pan_lin_overlap - shared_ids_before
    for mid in pan_lin_overlap:
        if mid in pan_map:
            del pan_map[mid]
        if mid in lin_map:
            del lin_map[mid]

    for mid in pan_lin_only_overlap:
        # build shared-style tuple from both callsets
        r_pan = pan_by_id.get(mid)
        r_lin = None
        # build lookup lazily from lin_rows_shared isn't enough; use full lin_effects
        # make a full exact lookup once
        if r_pan is None:
            continue
        # find linear exact record
        found = False
        chrom, pos = parse_marker_id(mid)
        for r in lin_effects.get((chrom, pos), []):
            if r["MarkerID"] == mid:
                r_lin = r
                found = True
                break
        if found and r_lin is not None:
            shared_final[mid] = (r_lin, r_pan)

    # --------------------------------
    # After-filter stats
    # --------------------------------
    final_pan_ids = set(pan_map)
    final_lin_ids = set(lin_map)
    final_shared_ids = set(shared_final)

    # final sanity
    overlap_pan_lin = final_pan_ids & final_lin_ids
    overlap_pan_shared = final_pan_ids & final_shared_ids
    overlap_lin_shared = final_lin_ids & final_shared_ids

    if overlap_pan_lin or overlap_pan_shared or overlap_lin_shared:
        sys.exit(
            "ERROR: Final sets still overlap: "
            f"pan∩lin={len(overlap_pan_lin)}, "
            f"pan∩shared={len(overlap_pan_shared)}, "
            f"lin∩shared={len(overlap_lin_shared)}"
        )

    print("After overlap filtering:")
    print(f"  Final Pangenie-unique IDs: {len(final_pan_ids):,}")
    print(f"  Final Linear-unique IDs:   {len(final_lin_ids):,}")
    print(f"  Final Shared pre-conc IDs: {len(final_shared_ids):,}")

    print("Moved during filtering:")
    print(f"  Removed from pan because also shared: {len(pan_shared_overlap):,}")
    print(f"  Removed from lin because also shared: {len(lin_shared_overlap):,}")
    print(f"  Removed from both unique sets and moved to shared: {len(pan_lin_only_overlap):,}")

    # --------------------------------
    # E) Write final outputs
    # --------------------------------
    pan_lines_final = list(pan_map.values())
    lin_lines_final = list(lin_map.values())

    with open(f"{prefix}.only_in_pangenie.tsv", "w") as f:
        f.write(header1 + "\n")
        for line in pan_lines_final:
            f.write(line + "\n")
    print(f"Wrote: {prefix}.only_in_pangenie.tsv")

    with open(f"{prefix}.only_in_linear.tsv", "w") as f:
        f.write(header2 + "\n")
        for line in lin_lines_final:
            f.write(line + "\n")
    print(f"Wrote: {prefix}.only_in_linear.tsv")

    with open(f"{prefix}.only_in_pangenie.ids.txt", "w") as f:
        for mid in sorted(final_pan_ids):
            f.write(mid + "\n")
    print(f"Wrote: {prefix}.only_in_pangenie.ids.txt")

    with open(f"{prefix}.only_in_linear.ids.txt", "w") as f:
        for mid in sorted(final_lin_ids):
            f.write(mid + "\n")
    print(f"Wrote: {prefix}.only_in_linear.ids.txt")

    # Write shared outputs in the same style as the other script
    write_shared_preconcordance_outputs(shared_final, prefix)

    # --------------------------------
    # Other summary stats
    # --------------------------------
    print(f"Linear-unique hap.py events:   {len(lin_pos):,}")
    print(f"Shared-pre hap.py events:      {len(shared_pos):,}")
    print(f"Pangenie-unique hap.py events: {len(pan_pos):,}")

    pan_happy_including_sv = len(pan_pos) + len(pan_skipped_sv_ids)
    print(f"Pangenie-unique hap.py events + skipped SV IDs: {pan_happy_including_sv:,}")
    print(f"Difference (hap.py+SV minus final recovered):   {pan_happy_including_sv - len(final_pan_ids):,}")

    print(f"Pangenie removed multisite IDs excluded: {len(pan_removed_multisite_ids):,}")
    print(f"Linear removed multisite IDs excluded:   {len(lin_removed_multisite_ids):,}")
    print(f"Pangenie skipped SV IDs requested for re-add: {len(pan_skipped_sv_ids):,}")
    print(f"Pangenie skipped SV IDs missing from GWAS:    {pan_sv_missing:,}")


    # Keep the old AF unique plots on the final unique sets
    plot_af_histograms(
        header1, idx1, pan_lines_final,
        header2, idx2, lin_lines_final,
        base1, base2, prefix
    )


if __name__ == "__main__":
    main()


