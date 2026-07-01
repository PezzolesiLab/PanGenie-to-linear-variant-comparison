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


MIN_P = 1e-300


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
    """
    Parse chrom + position from MarkerID.
    Works for:
      chrom_pos_ref_alt
      chrom_pos_ref_alt_zN
      chrom_pos_ins_size
      chrom_pos_del_size
      chrom_pos_complex_size
      ... as long as the first two fields are chrom and pos.
    """
    parts = marker_id.split("_")
    if len(parts) < 2:
        raise ValueError(f"Cannot parse MarkerID: {marker_id}")
    chrom = canonical_chrom(parts[0])
    pos1 = int(parts[1])
    return chrom, pos1


def parse_marker_id_ref_alt(marker_id: str):
    """
    Best-effort parse of allele-based IDs.

    Accepted:
      chrom_pos_ref_alt
      chrom_pos_ref_alt_zN

    Returns:
      (chrom, pos1, ref, alt) or None
    """
    parts = marker_id.split("_")
    if len(parts) == 4:
        try:
            chrom = canonical_chrom(parts[0])
            pos1 = int(parts[1])
        except ValueError:
            return None
        return chrom, pos1, parts[2], parts[3]

    if len(parts) == 5 and parts[4].startswith("z") and parts[4][1:].isdigit():
        try:
            chrom = canonical_chrom(parts[0])
            pos1 = int(parts[1])
        except ValueError:
            return None
        return chrom, pos1, parts[2], parts[3]

    return None


def safe_logp(p):
    if p <= 0:
        p = MIN_P
    elif p > 1:
        p = 1.0
    return -math.log10(p)


# -------------------------
# Load ID / p-value inputs
# -------------------------
def load_id_list(path):
    """
    Load MarkerIDs from either:
      - plain text (one ID per line)
      - TSV with MarkerID column
    """
    ids = set()
    with open_maybe_gzip(path, "rt") as f:
        first = f.readline().strip()
        if not first:
            return ids

        fields = first.split()

        if "MarkerID" in fields:
            mid_idx = fields.index("MarkerID")
            for line in f:
                if not line.strip():
                    continue
                parts = line.split()
                if len(parts) > mid_idx:
                    ids.add(parts[mid_idx])
        else:
            ids.add(first)
            for line in f:
                line = line.strip()
                if line:
                    ids.add(line)

    return ids


def load_gwas_pvalues(path):
    """
    Returns dict:
      MarkerID -> p.value
    """
    pvals = {}
    with open_maybe_gzip(path, "rt") as f:
        header = f.readline().split()
        try:
            mid_idx = header.index("MarkerID")
            p_idx = header.index("p.value")
        except ValueError as e:
            sys.exit(f"ERROR: required column missing in {path}: {e}")

        for line in f:
            if not line.strip():
                continue
            fields = line.split()
            try:
                pvals[fields[mid_idx]] = float(fields[p_idx])
            except Exception:
                continue

    return pvals


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
        total += bisect_left(sorted_positions, e) - bisect_left(sorted_positions, s)
    return total


# -------------------------
# Statistics
# -------------------------
def clamp_p(p):
    if not math.isfinite(p):
        return float("nan")
    return max(MIN_P, min(1.0, p))


def hypergeom_two_sided(N, K, n, k):
    if any(x < 0 for x in (N, K, n, k)) or K > N or k > n:
        return float("nan"), float("nan")

    if N <= 1:
        return 0.0, 1.0

    p0 = K / N
    mu = n * p0
    var = n * p0 * (1 - p0) * (N - n) / (N - 1)

    if var <= 0:
        return 0.0, 1.0

    z = (k - mu) / math.sqrt(var)
    p = 2 * min(
        0.5 * math.erfc(z / math.sqrt(2)),
        0.5 * math.erfc(-z / math.sqrt(2))
    )
    return z, clamp_p(p)


# -------------------------
# Concordance
# -------------------------
def load_concordance_table(path):
    """
    Accepts a concordance table with columns:
      Chromosome Start End Ref Alt Concordance(%)
    Returns:
      (chrom, pos1, ref, alt) -> concordance_pct
    """
    m = {}
    with open_maybe_gzip(path, "rt") as f:
        header = f.readline().split()
        idx = {c: header.index(c) for c in header}
        for line in f:
            if not line.strip():
                continue
            fields = line.split()
            try:
                chrom = canonical_chrom(fields[idx["Chromosome"]])
                pos1 = int(fields[idx["Start"]]) + 1
                ref = fields[idx["Ref"]]
                alt = fields[idx["Alt"]]
                conc = float(fields[idx["Concordance(%)"]])
            except Exception:
                continue
            m[(chrom, pos1, ref, alt)] = conc
    return m


def plot_concordance_hist(prefix, values):
    if not values:
        print("No concordance values to plot.")
        return
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(values, bins=30, range=(0, 100))
    ax.set_xlabel("Concordance (%)")
    ax.set_ylabel("Variant count")
    ax.set_title("Concordance of shared p-fold outliers")
    fig.tight_layout()
    out = f"{prefix}.concordance_hist.pdf"
    ensure_parent_dir(out)
    fig.savefig(out, format="pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote concordance histogram: {out}")


# -------------------------
# Plot enrichment
# -------------------------
def plot_enrichment(prefix, results, alpha, n, N):
    plot_rows = [r for r in results if math.isfinite(r["p_adj"]) and math.isfinite(r["fold"])]
    if not plot_rows:
        print("No valid rows to plot; skipping enrichment plot.")
        return

    plot_rows.sort(key=lambda r: (r["p_adj"] >= alpha, -r["fold"]))

    names = [r["bed"] for r in plot_rows]
    folds = [r["fold"] for r in plot_rows]
    fracs = [r["frac"] for r in plot_rows]
    p_adj = [r["p_adj"] for r in plot_rows]

    x = range(len(names))
    fig, ax = plt.subplots(figsize=(max(8, 0.45 * len(names)), 5))
    ax.bar(x, [f - 1 for f in folds], bottom=1.0)

    ax.axhline(1.0)
    ax.set_xticks(list(x))
    ax.set_xticklabels(names, rotation=90, fontsize=6)
    ax.set_ylabel("Fold enrichment (shared p-fold outliers vs analysis universe)")
    ax.set_title(
        f"Shared variants with large p-value changes\n"
        f"(background = pangenie-unique ∪ linear-unique ∪ shared; "
        f"excluded IDs removed; n={n:,}, N={N:,})"
    )

    yspan = ax.get_ylim()[1] - ax.get_ylim()[0]
    for i, (pa, fr, fval) in enumerate(zip(p_adj, fracs, folds)):
        if pa < alpha:
            top = max(fval, 1.0) if math.isfinite(fval) else 1.0
            ax.text(i, top + 0.05 * yspan, "*", ha="center", fontsize=11)
            ax.text(i, top + 0.025 * yspan, f"{fr*100:.1f}%", ha="center", fontsize=8)

    fig.tight_layout()
    out = f"{prefix}.shared_pfold_enrichment.pdf"
    ensure_parent_dir(out)
    fig.savefig(out, format="pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote enrichment plot: {out}")


# -------------------------
# Position conversion
# -------------------------
def ids_to_positions(ids):
    """
    Convert MarkerIDs to 0-based BED-like positions by chromosome.

    Returns:
      pos_by_chrom: dict[chrom] = sorted list of 0-based positions
      n_parsed
      n_failed
    """
    pos_by_chrom = defaultdict(list)
    n_parsed = 0
    n_failed = 0

    for mid in ids:
        try:
            ck, pos1 = parse_marker_id(mid)
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
            "Shared p-fold enrichment vs analysis-universe background.\n"
            "Background = pangenie-unique ∪ linear-unique ∪ shared.\n"
            "Excluded multisite IDs are removed from both background and outlier target."
        )
    )

    # analysis-universe set inputs
    ap.add_argument("--pangenie-unique", required=True,
                    help="Recovered pangenie-unique set (GWAS TSV with MarkerID column or plain ID list)")
    ap.add_argument("--linear-unique", required=True,
                    help="Recovered linear-unique set (GWAS TSV with MarkerID column or plain ID list)")
    ap.add_argument("--shared-ids", required=True,
                    help="Shared set used for this analysis (typically post-concordance)")

    # p-values still come from GWAS result tables
    ap.add_argument("pangenie_gwas")
    ap.add_argument("linear_gwas")

    # exclusions
    ap.add_argument("--excluded-ids", nargs="+", default=[],
                    help="ID lists to exclude from both background and outlier target, e.g. removed multisite IDs")

    # enrichment / plot controls
    ap.add_argument("--beds", nargs="+", required=True)
    ap.add_argument("--p-fold", type=float, default=10.0)
    ap.add_argument("--direction", choices=["pangenie", "linear", "both"], default="both")
    ap.add_argument("--concordance", default=None)
    ap.add_argument("--prefix", default="shared_pfold")
    ap.add_argument("--alpha", type=float, default=0.05)

    args = ap.parse_args()
    prefix = args.prefix

    # -------------------------
    # Load core sets
    # -------------------------
    ids_pan_unique_raw = load_id_list(args.pangenie_unique)
    ids_lin_unique_raw = load_id_list(args.linear_unique)
    ids_shared_raw = load_id_list(args.shared_ids)

    excluded_ids = set()
    for path in args.excluded_ids:
        excluded_ids |= load_id_list(path)

    # Background universe = actual analysis universe, not full GWAS union
    background_ids_raw = ids_pan_unique_raw | ids_lin_unique_raw | ids_shared_raw
    background_ids = background_ids_raw - excluded_ids

    # Shared set used for outlier selection
    shared_ids = ids_shared_raw - excluded_ids

    # Load p-values from full GWAS tables, but restrict by shared set
    p_pan = load_gwas_pvalues(args.pangenie_gwas)
    p_lin = load_gwas_pvalues(args.linear_gwas)

    # Shared IDs must have p-values in both callsets
    shared_ids &= p_pan.keys() & p_lin.keys()

    print(f"Pangenie-unique raw IDs:               {len(ids_pan_unique_raw):,}")
    print(f"Linear-unique raw IDs:                 {len(ids_lin_unique_raw):,}")
    print(f"Shared raw IDs:                        {len(ids_shared_raw):,}")
    print(f"Excluded IDs removed from analysis:    {len(excluded_ids):,}")
    print(f"Background union IDs after exclusion:  {len(background_ids):,}")
    print(f"Shared IDs after exclusion/pval check: {len(shared_ids):,}")

    if not background_ids:
        sys.exit("ERROR: Background universe is empty after exclusions.")
    if not shared_ids:
        sys.exit("ERROR: Shared set is empty after exclusions / p-value intersection.")

    # -------------------------
    # p-fold outlier selection inside the shared set
    # -------------------------
    thr = math.log10(args.p_fold)
    outlier_ids = set()

    for mid in shared_ids:
        dp = safe_logp(p_pan[mid]) - safe_logp(p_lin[mid])

        if args.direction == "pangenie" and dp >= thr:
            outlier_ids.add(mid)
        elif args.direction == "linear" and dp <= -thr:
            outlier_ids.add(mid)
        elif args.direction == "both" and abs(dp) >= thr:
            outlier_ids.add(mid)

    print(f"Shared p-fold outliers:                {len(outlier_ids):,}")

    if not outlier_ids:
        print("WARNING: No shared p-fold outliers found. Continuing to write empty/degenerate outputs where possible.")

    # -------------------------
    # Optional concordance histogram for outliers
    # -------------------------
    conc_vals = []
    if args.concordance:
        conc_map = load_concordance_table(args.concordance)
        for mid in outlier_ids:
            key = parse_marker_id_ref_alt(mid)
            if key and key in conc_map:
                conc_vals.append(conc_map[key])

        if conc_vals:
            print(f"Mean concordance of matched outliers:  {sum(conc_vals)/len(conc_vals):.2f}%")
        else:
            print("No concordance values matched outlier IDs.")

        plot_concordance_hist(prefix, conc_vals)

    # -------------------------
    # Build positions using analysis-universe background
    # -------------------------
    all_pos, N, N_failed = ids_to_positions(background_ids)
    tgt_pos, n, n_failed = ids_to_positions(outlier_ids)

    print(f"Background parsed positions (N):       {N:,}")
    print(f"Background IDs failed to parse:        {N_failed:,}")
    print(f"Outlier parsed positions (n):          {n:,}")
    print(f"Outlier IDs failed to parse:           {n_failed:,}")

    if N == 0:
        sys.exit("ERROR: No background IDs could be parsed into positions.")

    results = []
    for bed in args.beds:
        merged, nraw, nmerged, total_bases = load_bed_merged(bed)

        K = 0
        k = 0
        for ck, ivs in merged.items():
            K += count_positions_in_intervals(all_pos.get(ck, []), ivs)
            k += count_positions_in_intervals(tgt_pos.get(ck, []), ivs)

        frac = (k / n) if n else 0.0
        frac_bg = (K / N) if N else 0.0
        fold = (frac / frac_bg) if K and n and N else float("nan")
        z, p = hypergeom_two_sided(N, K, n, k) if n else (float("nan"), float("nan"))
        p_adj = clamp_p(p * len(args.beds)) if math.isfinite(p) else float("nan")

        results.append(
            {
                "bed": Path(bed).name,
                "bed_path": str(bed),
                "raw_intervals": nraw,
                "merged_intervals": nmerged,
                "covered_bases": total_bases,
                "k": k,
                "K": K,
                "frac": frac,
                "frac_bg": frac_bg,
                "fold": fold,
                "z": z,
                "p_adj": p_adj,
            }
        )

        print(
            f"{Path(bed).name}: "
            f"k={k:,}, K={K:,}, frac={frac:.6f}, frac_bg={frac_bg:.6f}, "
            f"fold={fold:.3f}, p_adj={p_adj:.3e}"
        )

    # -------------------------
    # TSV
    # -------------------------
    out_tsv = f"{prefix}.shared_pfold_enrichment.tsv"
    ensure_parent_dir(out_tsv)
    with open(out_tsv, "w") as f:
        f.write(
            "BedFile\tBedPath\tRawIntervals\tMergedIntervals\tCoveredBases\t"
            "N_background\tn_outliers\tK_background_in_bed\tk_outliers_in_bed\t"
            "Frac_background\tFrac_outliers\tFold\tZ\tP_adj\n"
        )
        for r in results:
            f.write(
                f"{r['bed']}\t{r['bed_path']}\t{r['raw_intervals']}\t{r['merged_intervals']}\t{r['covered_bases']}\t"
                f"{N}\t{n}\t{r['K']}\t{r['k']}\t"
                f"{r['frac_bg']:.6f}\t{r['frac']:.6f}\t{r['fold']:.6f}\t{r['z']:.6f}\t{r['p_adj']:.3e}\n"
            )
    print(f"Wrote TSV: {out_tsv}")

    # -------------------------
    # Plot
    # -------------------------
    plot_enrichment(prefix, results, args.alpha, n, N)


if __name__ == "__main__":
    main()







