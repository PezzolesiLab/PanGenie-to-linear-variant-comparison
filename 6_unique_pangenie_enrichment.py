#!/usr/bin/env python3
import argparse
import gzip
import sys
from pathlib import Path
import math
from bisect import bisect_left
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


MIN_P = 1e-300  # floor for extremely small p-values


# -------------------------
# IO helpers
# -------------------------
def open_maybe_gzip(path, mode="rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def ensure_parent_dir(pathlike):
    p = Path(pathlike)
    if p.parent and str(p.parent) not in (".", ""):
        p.parent.mkdir(parents=True, exist_ok=True)


def canonical_chrom(chrom: str) -> str:
    chrom = chrom.strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    return chrom.lower()


def parse_marker_id(marker_id: str):
    parts = marker_id.split("_")
    if len(parts) < 2:
        raise ValueError(f"Cannot parse MarkerID: {marker_id}")
    chrom_raw = parts[0]
    pos1 = int(parts[1])
    return chrom_raw, pos1, canonical_chrom(chrom_raw)


# -------------------------
# Load IDs
# -------------------------
def load_id_list(path):
    """
    Load MarkerIDs from either:
      - a plain text file (one MarkerID per line), OR
      - a GWAS-style TSV containing a 'MarkerID' column.

    Returns: set of MarkerID strings.
    """
    ids = set()

    with open_maybe_gzip(path, "rt") as f:
        first = f.readline().rstrip("\n")
        if not first:
            return ids

        fields = first.split()

        # Case 1: TSV with MarkerID column
        if "MarkerID" in fields:
            mid_idx = fields.index("MarkerID")
            for line in f:
                if not line.strip():
                    continue
                parts = line.split()
                if len(parts) <= mid_idx:
                    continue
                ids.add(parts[mid_idx])

        # Case 2: plain ID list
        else:
            ids.add(first.strip())
            for line in f:
                line = line.strip()
                if line:
                    ids.add(line)

    return ids


# -------------------------
# BED utilities
# -------------------------
def merge_intervals(intervals):
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = []
    s0, e0 = intervals[0]
    for s, e in intervals[1:]:
        if s <= e0:
            e0 = max(e0, e)
        else:
            merged.append((s0, e0))
            s0, e0 = s, e
    merged.append((s0, e0))
    return merged


def load_bed_merged(path):
    raw = defaultdict(list)
    num_raw = 0

    with open_maybe_gzip(path, "rt") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, s, e = line.split()[:3]
            ck = canonical_chrom(chrom)
            try:
                s = int(s)
                e = int(e)
            except ValueError:
                continue
            raw[ck].append((s, e))
            num_raw += 1

    merged = {}
    total_bases = 0
    for ck, ivs in raw.items():
        m = merge_intervals(ivs)
        merged[ck] = m
        for s, e in m:
            total_bases += (e - s)

    return merged, num_raw, sum(len(v) for v in merged.values()), total_bases


def count_positions_in_intervals(sorted_positions, intervals):
    total = 0
    for s, e in intervals:
        left = bisect_left(sorted_positions, s)
        right = bisect_left(sorted_positions, e)
        total += (right - left)
    return total


# -------------------------
# Statistics
# -------------------------
def clamp_p(p):
    if not math.isfinite(p):
        return float("nan")
    return max(MIN_P, min(1.0, p))


def hypergeom_normal_pvalues_two_sided(N, K, n, k):
    if any(x < 0 for x in (N, K, n, k)) or K > N or k > n:
        return float("nan"), float("nan"), float("nan"), float("nan")

    p0 = K / N
    mu = n * p0
    var = n * p0 * (1 - p0) * (N - n) / (N - 1) if N > 1 else 0.0

    if var <= 0:
        return 0.0, 1.0, 1.0, 1.0

    z = (k - mu) / math.sqrt(var)
    p_right = 0.5 * math.erfc(z / math.sqrt(2))
    p_left = 0.5 * math.erfc(-z / math.sqrt(2))
    p_two = clamp_p(2 * min(p_right, p_left))
    return z, clamp_p(p_right), clamp_p(p_left), p_two


# -------------------------
# Position conversion
# -------------------------
def ids_to_positions(ids):
    """
    Convert MarkerIDs to 0-based BED-like positions by chromosome.

    Returns:
      pos_by_chrom: dict[chrom_key] = sorted list of 0-based positions
      n_parsed: number of IDs successfully parsed
      n_failed: number of IDs that failed parsing
    """
    pos_by_chrom = defaultdict(list)
    n_parsed = 0
    n_failed = 0

    for mid in ids:
        try:
            _, pos1, ck = parse_marker_id(mid)
            pos_by_chrom[ck].append(pos1 - 1)
            n_parsed += 1
        except Exception:
            n_failed += 1

    for ck in pos_by_chrom:
        pos_by_chrom[ck].sort()

    return pos_by_chrom, n_parsed, n_failed


# -------------------------
# Main
# -------------------------
def main():
    ap = argparse.ArgumentParser(
        description=(
            "Pangenie-unique enrichment using the analysis universe defined by:\n"
            "  pangenie-unique ∪ linear-unique ∪ shared (pre-concordance),\n"
            "with removed multisite IDs excluded from both background and target."
        )
    )

    # Core set inputs
    ap.add_argument("--pangenie-unique", required=True,
                    help="Recovered pangenie-unique set (GWAS TSV with MarkerID column or plain ID list)")
    ap.add_argument("--linear-unique", required=True,
                    help="Recovered linear-unique set (GWAS TSV with MarkerID column or plain ID list)")
    ap.add_argument("--shared-pre", required=True,
                    help="Shared pre-concordance set (GWAS TSV with MarkerID column or plain ID list)")

    # Exclusions
    ap.add_argument("--excluded-ids", nargs="+", default=[],
                    help="ID lists to exclude from both background and target, e.g. removed multisite IDs")

    # Region tests
    ap.add_argument("--beds", nargs="+", required=True)
    ap.add_argument("--prefix", default="pangenie_unique")
    ap.add_argument("--alpha", type=float, default=0.05)

    args = ap.parse_args()

    prefix = args.prefix

    # -------------------------
    # Load sets
    # -------------------------
    ids_pan_unique_raw = load_id_list(args.pangenie_unique)
    ids_lin_unique_raw = load_id_list(args.linear_unique)
    ids_shared_pre_raw = load_id_list(args.shared_pre)

    excluded_ids = set()
    for path in args.excluded_ids:
        excluded_ids |= load_id_list(path)

    # Background universe = union of recovered unique/shared sets, excluding multisite IDs
    all_ids_raw = ids_pan_unique_raw | ids_lin_unique_raw | ids_shared_pre_raw
    all_ids = all_ids_raw - excluded_ids

    # Target = pangenie-unique only, excluding multisite IDs
    unique_ids = ids_pan_unique_raw - excluded_ids

    print(f"Pangenie-unique raw IDs:                 {len(ids_pan_unique_raw):,}")
    print(f"Linear-unique raw IDs:                   {len(ids_lin_unique_raw):,}")
    print(f"Shared pre-concordance raw IDs:          {len(ids_shared_pre_raw):,}")
    print(f"Excluded IDs removed from analysis:      {len(excluded_ids):,}")
    print(f"Background union IDs after exclusion:    {len(all_ids):,}")
    print(f"Pangenie-unique IDs after exclusion:     {len(unique_ids):,}")

    if not all_ids:
        sys.exit("ERROR: Background union is empty after exclusions.")
    if not unique_ids:
        sys.exit("ERROR: Pangenie-unique target set is empty after exclusions.")

    # Convert to positions
    all_pos, N, N_failed = ids_to_positions(all_ids)
    unique_pos, n, n_failed = ids_to_positions(unique_ids)

    print(f"Background parsed positions (N):         {N:,}")
    print(f"Background IDs failed to parse:          {N_failed:,}")
    print(f"Unique parsed positions (n):             {n:,}")
    print(f"Unique IDs failed to parse:              {n_failed:,}")

    if N == 0:
        sys.exit("ERROR: No background IDs could be parsed into positions.")
    if n == 0:
        sys.exit("ERROR: No unique IDs could be parsed into positions.")

    # -------------------------
    # Per-BED enrichment
    # -------------------------
    results = []
    m_tests = len(args.beds)

    for bed in args.beds:
        merged, nraw, nmerged, total_bases = load_bed_merged(bed)

        K = 0
        k = 0
        for ck, ivs in merged.items():
            K += count_positions_in_intervals(all_pos.get(ck, []), ivs)
            k += count_positions_in_intervals(unique_pos.get(ck, []), ivs)

        z, p_r, p_l, p_two = hypergeom_normal_pvalues_two_sided(N, K, n, k)
        p_adj = clamp_p(p_two * m_tests)

        frac_bg = (K / N) if N > 0 else float("nan")
        frac_unique = (k / n) if n > 0 else float("nan")
        fold = (frac_unique / frac_bg) if K > 0 and N > 0 else float("nan")

        results.append(
            {
                "bed": Path(bed).name,
                "bed_path": str(bed),
                "nraw_intervals": nraw,
                "nmerged_intervals": nmerged,
                "total_bases": total_bases,
                "K": K,
                "k": k,
                "frac_bg": frac_bg,
                "frac_unique": frac_unique,
                "fold": fold,
                "z": z,
                "p_right": p_r,
                "p_left": p_l,
                "p_two": p_two,
                "p_adj": p_adj,
            }
        )

        print(
            f"{Path(bed).name}: "
            f"k={k:,}, K={K:,}, frac_unique={frac_unique:.6f}, frac_bg={frac_bg:.6f}, "
            f"fold={fold:.3f}, p_two={p_two:.3e}, p_adj={p_adj:.3e}"
        )

    # -------------------------
    # Write TSV
    # -------------------------
    out_tsv = f"{prefix}.pangenie_unique_enrichment.tsv"
    ensure_parent_dir(out_tsv)
    with open(out_tsv, "w") as f:
        f.write(
            "BedFile\tBedPath\tRawIntervals\tMergedIntervals\tCoveredBases\t"
            "N_background\tn_unique\tK_background_in_bed\tk_unique_in_bed\t"
            "Frac_background\tFrac_unique\tFold\tZ\tP_right\tP_left\tP_two\tP_adj\n"
        )
        for r in results:
            f.write(
                f"{r['bed']}\t{r['bed_path']}\t{r['nraw_intervals']}\t{r['nmerged_intervals']}\t{r['total_bases']}\t"
                f"{N}\t{n}\t{r['K']}\t{r['k']}\t"
                f"{r['frac_bg']}\t{r['frac_unique']}\t{r['fold']}\t{r['z']}\t"
                f"{r['p_right']}\t{r['p_left']}\t{r['p_two']}\t{r['p_adj']}\n"
            )

    print(f"Wrote: {out_tsv}")

    # -------------------------
    # Plot fold enrichment (PDF output)
    # -------------------------
    plot_rows = [r for r in results if math.isfinite(r["p_adj"]) and math.isfinite(r["fold"])]
    if not plot_rows:
        print("No valid adjusted p-values / fold values to plot; skipping plot.")
        return

    def is_sig(r):
        return r["p_adj"] < args.alpha

    # Sort: significant first, then by fold (descending)
    def sort_key(r):
        sig_group = 0 if is_sig(r) else 1
        fold_val = r["fold"]
        if not math.isfinite(fold_val):
            fold_val = float("-inf")
        return (sig_group, -fold_val)

    plot_rows.sort(key=sort_key)

    names = [r["bed"] for r in plot_rows]
    folds = [r["fold"] for r in plot_rows]
    p_adj = [r["p_adj"] for r in plot_rows]
    frac_unique = [r["frac_unique"] for r in plot_rows]

    x = list(range(len(names)))
    fig, ax = plt.subplots(figsize=(max(8, 0.45 * len(names)), 5))

    heights = [(f - 1.0) if math.isfinite(f) else 0.0 for f in folds]
    bars = ax.bar(x, heights, bottom=1.0, alpha=0.85)

    ax.axhline(1.0, linewidth=1.0)

    ax.set_ylabel("Fold enrichment (Pangenie-unique vs analysis universe)", fontsize=10)
    ax.set_xlabel("Genomic region (BED)", fontsize=10)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=90, fontsize=6)
    ax.set_title(
        f"Pangenie-unique enrichment by genomic region\n"
        f"(background = pangenie-unique ∪ linear-unique ∪ shared-pre; "
        f"excluded IDs removed; two-sided test, Bonferroni; n={n:,}, N={N:,})",
        fontsize=12,
    )

    finite_folds = [f for f in folds if math.isfinite(f)]
    if finite_folds:
        ymin = min(finite_folds + [1.0])
        ymax = max(finite_folds + [1.0])
        pad = (ymax - ymin) * 0.15 if ymax > ymin else 0.2
        ax.set_ylim(ymin - pad, ymax + pad)

    yspan = ax.get_ylim()[1] - ax.get_ylim()[0]
    star_offset = 0.04 * yspan
    pct_offset = 0.015 * yspan

    for bar, pa, fu, fval in zip(bars, p_adj, frac_unique, folds):
        if pa >= args.alpha:
            continue

        bar_top = max(fval, 1.0) if math.isfinite(fval) else 1.0

        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            bar_top + star_offset,
            "*",
            ha="center",
            va="bottom",
            fontsize=11,
            fontweight="bold",
        )

        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            bar_top + pct_offset,
            f"{fu * 100:.1f}%",
            ha="center",
            va="bottom",
            fontsize=8,
        )

    fig.tight_layout()
    plot_path = f"{prefix}.pangenie_unique_twosided_fold_plot.pdf"
    ensure_parent_dir(plot_path)
    fig.savefig(plot_path, format="pdf", bbox_inches="tight")
    plt.close(fig)

    print(f"Fold enrichment plot written to: {plot_path}")


if __name__ == "__main__":
    main()






