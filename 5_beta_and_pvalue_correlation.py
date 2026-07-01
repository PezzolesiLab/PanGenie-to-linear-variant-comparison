#!/usr/bin/env python3
import argparse
import gzip
import sys
from pathlib import Path
from collections import defaultdict
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


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
    MarkerID must begin with chrom_pos_...
    Returns (chrom, pos).
    """
    parts = marker_id.split("_")
    if len(parts) < 2:
        raise ValueError
    chrom = canonical_chrom(parts[0])
    pos = int(parts[1])
    return chrom, pos


def parse_marker_id_with_alleles(marker_id):
    """
    Best-effort parse:
      - chrom_pos_ref_alt
      - chrom_pos_ref_alt_zN
    Returns (chrom, pos, ref, alt), or ref/alt=None if not parseable as allele-based ID.
    """
    parts = marker_id.split("_")
    if len(parts) < 2:
        raise ValueError

    chrom = canonical_chrom(parts[0])
    pos = int(parts[1])

    if len(parts) == 4:
        ref = parts[2]
        alt = parts[3]
        return chrom, pos, ref, alt

    if len(parts) == 5 and parts[4].startswith("z") and parts[4][1:].isdigit():
        ref = parts[2]
        alt = parts[3]
        return chrom, pos, ref, alt

    return chrom, pos, None, None


def safe_neg_log10(p):
    if p <= 0:
        p = 1e-300
    elif p > 1:
        p = 1.0
    return -math.log10(p)


def ensure_parent_dir(pathlike):
    p = Path(pathlike)
    if p.parent and str(p.parent) not in (".", ""):
        p.parent.mkdir(parents=True, exist_ok=True)


def safe_zscore(x):
    x = np.asarray(x, dtype=float)
    if x.size == 0:
        return x
    mu = np.nanmean(x)
    sd = np.nanstd(x, ddof=1)
    if not np.isfinite(sd) or sd == 0:
        return np.zeros_like(x)
    return (x - mu) / sd


def save_fig_both(fig, out_pdf, *, dpi=300):
    """
    Save the same figure to:
      - out_pdf (PDF)
      - out_pdf with .png extension (PNG)
    """
    ensure_parent_dir(out_pdf)
    out_pdf = str(out_pdf)
    if not out_pdf.lower().endswith(".pdf"):
        raise ValueError(f"Expected .pdf output path, got: {out_pdf}")

    out_png = out_pdf[:-4] + ".png"

    fig.savefig(out_pdf, format="pdf", bbox_inches="tight")
    fig.savefig(out_png, format="png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)


# -------------------------
# Load pre-concordance shared metrics
# -------------------------
def load_shared_pre_metrics(path):
    """
    Load the already-created pre-concordance shared metrics table.

    Expected columns:
      MarkerID
      BETA_pan SE_pan p_pan AF_pan
      BETA_lin SE_lin p_lin AF_lin

    Returns:
      header, records
    """
    path = Path(path)
    with open_maybe_gzip(path, "rt") as f:
        header = f.readline().rstrip("\n")
        if not header:
            sys.exit(f"ERROR: {path} empty")

        cols = header.split()
        req = [
            "MarkerID",
            "BETA_pan", "SE_pan", "p_pan", "AF_pan",
            "BETA_lin", "SE_lin", "p_lin", "AF_lin",
        ]
        try:
            idx = {c: cols.index(c) for c in req}
        except ValueError as e:
            sys.exit(f"ERROR: missing required column in {path}: {e}")

        records = []
        for line in f:
            if not line.strip():
                continue
            fields = line.split()
            try:
                rec = {
                    "MarkerID": fields[idx["MarkerID"]],
                    "BETA_pan": float(fields[idx["BETA_pan"]]),
                    "SE_pan": float(fields[idx["SE_pan"]]),
                    "p_pan": float(fields[idx["p_pan"]]),
                    "AF_pan": float(fields[idx["AF_pan"]]),
                    "BETA_lin": float(fields[idx["BETA_lin"]]),
                    "SE_lin": float(fields[idx["SE_lin"]]),
                    "p_lin": float(fields[idx["p_lin"]]),
                    "AF_lin": float(fields[idx["AF_lin"]]),
                }
            except Exception:
                continue
            records.append(rec)

    return header, records


# -------------------------
# Concordance BED loader
# -------------------------
def read_concordance_bed(path):
    """
    BED columns expected:
    Chromosome  Start  End  Ref  Alt  Concordance(%)  Total_Genotypes  Concordant_Genotypes

    Assumes BED Start is 0-based; converts to 1-based position using pos = Start + 1.
    Stores:
      - full key: (chrom,pos,ref,alt) -> (conc_pct,total,concordant)
      - pos key:  (chrom,pos) -> list of (ref,alt,conc_pct,total,concordant)
      - conc100_pos: set of (chrom,pos) where any record has conc==100
    """
    full = {}
    by_pos = defaultdict(list)
    conc100_pos = set()

    with open_maybe_gzip(path, "rt") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            low = line.lower()
            if low.startswith("chromosome") or low.startswith("#"):
                continue

            fields = line.split()
            if len(fields) < 8:
                continue

            chrom_raw, start_s, end_s, ref, alt, conc_s, total_s, concord_s = fields[:8]
            try:
                chrom = canonical_chrom(chrom_raw)
                start = int(start_s)
                pos = start + 1
                conc = float(conc_s)
                total = int(total_s)
                concord = int(concord_s)
            except Exception:
                continue

            full[(chrom, pos, ref, alt)] = (conc, total, concord)
            by_pos[(chrom, pos)].append((ref, alt, conc, total, concord))
            if abs(conc - 100.0) < 1e-9:
                conc100_pos.add((chrom, pos))

    return full, by_pos, conc100_pos


def lookup_concordance(marker_id, conc_full, conc_by_pos):
    """
    Returns (conc_pct,total,concordant) or None.

    Preference:
      1) If MarkerID has REF/ALT parseable, try full key (chrom,pos,ref,alt).
      2) Otherwise, if exactly one concordance record exists at (chrom,pos), use it.
      3) Otherwise None (ambiguous or missing).
    """
    chrom, pos, ref, alt = parse_marker_id_with_alleles(marker_id)
    if ref is not None and alt is not None:
        v = conc_full.get((chrom, pos, ref, alt))
        if v is not None:
            return v

    recs = conc_by_pos.get((chrom, pos), [])
    if len(recs) == 1:
        _, _, conc, total, concord = recs[0]
        return (conc, total, concord)

    return None


# -------------------------
# Plot helpers (PDF+PNG output)
# -------------------------
def scatter_plot(x, y, xlabel, ylabel, title, out_pdf, *, dpi=300):
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(x, y, s=5, alpha=0.5)
    mn = min(np.nanmin(x), np.nanmin(y))
    mx = max(np.nanmax(x), np.nanmax(y))
    ax.plot([mn, mx], [mn, mx], linestyle="--")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    fig.tight_layout()
    save_fig_both(fig, out_pdf, dpi=dpi)


def hexbin_plot(x, y, xlabel, ylabel, title, out_pdf, *, dpi=300, gridsize=60, bins="log"):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    fig, ax = plt.subplots(figsize=(6, 6))

    xy = np.concatenate([x.ravel(), y.ravel()])
    xy = xy[np.isfinite(xy)]
    if xy.size == 0:
        mn, mx = 0.0, 1.0
    else:
        mn, mx = float(np.min(xy)), float(np.max(xy))
        if mn == mx:
            mx = mn + 1.0

    hb = ax.hexbin(
        x, y,
        gridsize=gridsize,
        mincnt=1,
        bins=bins,
        cmap="Blues",
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    cbar = fig.colorbar(hb, ax=ax)
    cbar.set_label("log10(N)" if bins == "log" else "N")

    fig.tight_layout()
    save_fig_both(fig, out_pdf, dpi=dpi)


def hist_plot(vals, xlabel, title, out_pdf, *, dpi=300):
    fig, ax = plt.subplots(figsize=(8, 4))
    vmax = np.nanmax(np.abs(vals)) if len(vals) else 1.0
    if not np.isfinite(vmax) or vmax == 0:
        vmax = 1.0
    ax.hist(vals, bins=50, range=(-vmax, vmax))
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Count")
    ax.set_title(title)
    fig.tight_layout()
    save_fig_both(fig, out_pdf, dpi=dpi)


def z_plot(z, ylabel, title, out_pdf, *, dpi=300):
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.scatter(np.arange(len(z)), z, s=3, alpha=0.5)
    ax.axhline(0)
    for t in [2, 3, 4]:
        ax.axhline(t, linestyle="--")
        ax.axhline(-t, linestyle="--")
    ax.set_ylabel(ylabel)
    ax.set_xlabel("Shared variant index")
    ax.set_title(title)
    fig.tight_layout()
    save_fig_both(fig, out_pdf, dpi=dpi)


# -------------------------
# Main
# -------------------------
def main():
    ap = argparse.ArgumentParser(
        description=(
            "Use an existing pre-concordance shared metrics file, add concordance annotations, "
            "write post-concordance shared ID list, generate plots/metrics, and write a sidecar "
            "file for missing/ambiguous concordance variants."
        )
    )
    ap.add_argument(
        "--shared-pre-metrics",
        required=True,
        help="Existing pre-concordance shared metrics TSV"
    )
    ap.add_argument("--prefix", default="hap_shared")

    ap.add_argument("--concordance-bed", required=True,
                    help="BED with concordance percentages per variant")
    ap.add_argument("--min-concordance", type=float, default=100.0,
                    help="Keep variants with concordance >= this (default 100.0) for downstream shared ID list")

    ap.add_argument("--plot-dpi", type=int, default=300,
                    help="PNG DPI for plots (PDF is vector). Default 300.")

    args = ap.parse_args()

    # Load already-created pre-concordance shared metrics
    _, shared = load_shared_pre_metrics(args.shared_pre_metrics)
    print(f"Shared variants pre-concordance (loaded from file): {len(shared):,}")

    if not shared:
        sys.exit("No shared variants found in pre-concordance metrics file")

    # Load concordance and annotate/filter
    conc_full, conc_by_pos, conc100_pos = read_concordance_bed(args.concordance_bed)

    # preconcordance plot/metrics set = ALL shared variants from the provided file
    # each item:
    #   (rec, conc_pct_or_None, total_or_None, concord_or_None, passes_threshold_bool)
    pre_metrics = []

    # downstream filtered set = only post-concordance kept variants
    kept = []

    # sidecar rows for missing / ambiguous concordance
    missing_conc_rows = []

    missing_conc = 0
    below_thresh = 0

    for rec in shared:
        mid = rec["MarkerID"]
        chrom, pos, ref, alt = parse_marker_id_with_alleles(mid)
        recs_at_pos = conc_by_pos.get((chrom, pos), [])

        conc_rec = lookup_concordance(mid, conc_full, conc_by_pos)

        # optional position-only 100% fallback, same as previous behavior
        used_pos100_fallback = False
        if conc_rec is None:
            if (chrom, pos) in conc100_pos and args.min_concordance <= 100.0:
                conc_rec = (100.0, -1, -1)
                used_pos100_fallback = True

        if conc_rec is None:
            pre_metrics.append((rec, None, None, None, False))
            missing_conc += 1

            if recs_at_pos:
                bed_records = ";".join(
                    f"{rref},{ralt},{rconc},{rtotal},{rconcord}"
                    for (rref, ralt, rconc, rtotal, rconcord) in recs_at_pos
                )
            else:
                bed_records = "NA"

            missing_conc_rows.append({
                "MarkerID": mid,
                "Chrom": chrom,
                "Pos": pos,
                "Ref": ref if ref is not None else "NA",
                "Alt": alt if alt is not None else "NA",
                "ConcordanceRecordsAtPos": len(recs_at_pos),
                "UsedPos100Fallback": 1 if used_pos100_fallback else 0,
                "BedRecordsAtPos": bed_records,
                "BETA_pan": rec["BETA_pan"],
                "SE_pan": rec["SE_pan"],
                "p_pan": rec["p_pan"],
                "AF_pan": rec["AF_pan"],
                "BETA_lin": rec["BETA_lin"],
                "SE_lin": rec["SE_lin"],
                "p_lin": rec["p_lin"],
                "AF_lin": rec["AF_lin"],
            })
            continue

        conc_pct, total_gt, concord_gt = conc_rec
        passes = (conc_pct + 1e-9 >= args.min_concordance)
        pre_metrics.append((rec, conc_pct, total_gt, concord_gt, passes))

        if passes:
            kept.append((rec, conc_pct, total_gt, concord_gt))
        else:
            below_thresh += 1

    print("Concordance filtering stats:")
    print(f"  Shared variants pre-concordance (used for plots/metrics): {len(shared):,}")
    print(f"  Missing/ambiguous concordance:                            {missing_conc:,}")
    print(f"  Below threshold (<{args.min_concordance:g}%):             {below_thresh:,}")
    print(f"  Shared variants post-concordance (used downstream):       {len(kept):,}")

    # Write missing/ambiguous concordance variants sidecar
    missing_out = f"{args.prefix}.missing_or_ambiguous_concordance.tsv"
    ensure_parent_dir(missing_out)
    with open(missing_out, "w") as f:
        f.write(
            "MarkerID\tChrom\tPos\tRef\tAlt\tConcordanceRecordsAtPos\tUsedPos100Fallback\tBedRecordsAtPos\t"
            "BETA_pan\tSE_pan\tp_pan\tAF_pan\tBETA_lin\tSE_lin\tp_lin\tAF_lin\n"
        )
        for row in missing_conc_rows:
            f.write(
                f"{row['MarkerID']}\t{row['Chrom']}\t{row['Pos']}\t{row['Ref']}\t{row['Alt']}\t"
                f"{row['ConcordanceRecordsAtPos']}\t{row['UsedPos100Fallback']}\t{row['BedRecordsAtPos']}\t"
                f"{row['BETA_pan']}\t{row['SE_pan']}\t{row['p_pan']}\t{row['AF_pan']}\t"
                f"{row['BETA_lin']}\t{row['SE_lin']}\t{row['p_lin']}\t{row['AF_lin']}\n"
            )
    print(f"Wrote missing/ambiguous concordance variants: {missing_out}")

    if not kept:
        print("WARNING: No shared variants remain after concordance filtering for downstream analyses.")

    # Output ID list (POST-filter; main downstream list)
    shared_ids_path = f"{args.prefix}.shared_gwas_ids.txt"
    ensure_parent_dir(shared_ids_path)
    with open(shared_ids_path, "w") as f:
        for rec, _, _, _ in kept:
            f.write(rec["MarkerID"] + "\n")
    print(f"Wrote concordance-filtered shared ID list: {shared_ids_path}")

    # Arrays for plots/metrics come from the FULL pre-concordance shared set
    beta_lin = np.array([r["BETA_lin"] for r, _, _, _, _ in pre_metrics], dtype=float)
    beta_pan = np.array([r["BETA_pan"] for r, _, _, _, _ in pre_metrics], dtype=float)
    se_lin = np.array([r["SE_lin"] for r, _, _, _, _ in pre_metrics], dtype=float)
    se_pan = np.array([r["SE_pan"] for r, _, _, _, _ in pre_metrics], dtype=float)
    logp_lin = np.array([safe_neg_log10(r["p_lin"]) for r, _, _, _, _ in pre_metrics], dtype=float)
    logp_pan = np.array([safe_neg_log10(r["p_pan"]) for r, _, _, _, _ in pre_metrics], dtype=float)
    af_lin = np.array([r["AF_lin"] for r, _, _, _, _ in pre_metrics], dtype=float)
    af_pan = np.array([r["AF_pan"] for r, _, _, _, _ in pre_metrics], dtype=float)

    # Correlations on pre-concordance shared set
    if len(pre_metrics) >= 2:
        print(f"Pearson r BETA (pre-concordance): {np.corrcoef(beta_lin, beta_pan)[0,1]:.4f}")
        print(f"Pearson r logP (pre-concordance): {np.corrcoef(logp_lin, logp_pan)[0,1]:.4f}")
        print(f"Pearson r AF (pre-concordance):   {np.corrcoef(af_lin, af_pan)[0,1]:.4f}")
    else:
        print("Pearson r: n<2 in pre-concordance shared set; correlations not meaningful.")

    # Deltas + z on pre-concordance shared set
    delta_beta = beta_pan - beta_lin
    z_beta = delta_beta / np.sqrt(se_pan**2 + se_lin**2)

    delta_logp = logp_pan - logp_lin
    z_logp = safe_zscore(delta_logp)

    delta_af = af_pan - af_lin
    z_af = safe_zscore(delta_af)

    dpi = args.plot_dpi

    # Plots (PDF+PNG) from pre-concordance shared set
    scatter_plot(
        beta_pan, beta_lin,
        "BETA (Pangenie)", "BETA (Linear)",
        "Shared variants (pre-concordance): BETA comparison",
        f"{args.prefix}.shared_beta_scatter.pdf",
        dpi=dpi
    )

    scatter_plot(
        logp_pan, logp_lin,
        "-log10(p) (Pangenie)", "-log10(p) (Linear)",
        "Shared variants (pre-concordance): -log10(p) comparison",
        f"{args.prefix}.shared_logp_scatter.pdf",
        dpi=dpi
    )

    scatter_plot(
        af_pan, af_lin,
        "AF (Pangenie)", "AF (Linear)",
        "Shared variants (pre-concordance): AF comparison",
        f"{args.prefix}.shared_af_scatter.pdf",
        dpi=dpi
    )

    hexbin_plot(
        beta_pan, beta_lin,
        "BETA (Pangenie)", "BETA (Linear)",
        "Shared variants (pre-concordance): BETA hexbin density",
        f"{args.prefix}.shared_beta_hexbin.pdf",
        dpi=dpi
    )

    hexbin_plot(
        logp_pan, logp_lin,
        "-log10(p) (Pangenie)", "-log10(p) (Linear)",
        "Shared variants (pre-concordance): -log10(p) hexbin density",
        f"{args.prefix}.shared_logp_hexbin.pdf",
        dpi=dpi
    )

    hexbin_plot(
        af_pan, af_lin,
        "AF (Pangenie)", "AF (Linear)",
        "Shared variants (pre-concordance): AF hexbin density",
        f"{args.prefix}.shared_af_hexbin.pdf",
        dpi=dpi
    )

    hist_plot(
        delta_beta, "ΔBETA (Pangenie − Linear)",
        "Shared variants (pre-concordance): ΔBETA",
        f"{args.prefix}.shared_delta_beta_hist.pdf",
        dpi=dpi
    )

    hist_plot(
        delta_logp, "Δ−log10(p)",
        "Shared variants (pre-concordance): Δ−log10(p)",
        f"{args.prefix}.shared_delta_logp_hist.pdf",
        dpi=dpi
    )

    hist_plot(
        delta_af, "ΔAF",
        "Shared variants (pre-concordance): ΔAF",
        f"{args.prefix}.shared_delta_af_hist.pdf",
        dpi=dpi
    )

    z_plot(
        z_beta, "z_diff(beta)",
        "Shared variants (pre-concordance): z_diff(beta)",
        f"{args.prefix}.shared_z_beta.pdf",
        dpi=dpi
    )

    z_plot(
        z_logp, "z_diff(logp)",
        "Shared variants (pre-concordance): z_diff(logp)",
        f"{args.prefix}.shared_z_logp.pdf",
        dpi=dpi
    )

    z_plot(
        z_af, "z_diff(AF)",
        "Shared variants (pre-concordance): z_diff(AF)",
        f"{args.prefix}.shared_z_af.pdf",
        dpi=dpi
    )

    # Metrics table on pre-concordance shared set, now including concordance annotations
    out_tsv = f"{args.prefix}.shared_variants_metrics.pre_concordance.tsv"
    ensure_parent_dir(out_tsv)
    with open(out_tsv, "w") as f:
        f.write(
            "MarkerID\t"
            "Concordance_pct\tTotal_Genotypes\tConcordant_Genotypes\tPassesConcordanceThreshold\t"
            "BETA_pan\tSE_pan\tp_pan\tAF_pan\t"
            "BETA_lin\tSE_lin\tp_lin\tAF_lin\t"
            "delta_beta\tz_beta\tdelta_logp\tz_logp\tdelta_af\tz_af\n"
        )
        for (r, conc_pct, total_gt, concord_gt, passes), db, zb, dlp, zlp, daf, zaf in zip(
            pre_metrics, delta_beta, z_beta, delta_logp, z_logp, delta_af, z_af
        ):
            conc_pct_str = "NA" if conc_pct is None else str(conc_pct)
            total_gt_str = "NA" if total_gt is None else str(total_gt)
            concord_gt_str = "NA" if concord_gt is None else str(concord_gt)
            passes_str = "1" if passes else "0"

            f.write(
                f"{r['MarkerID']}\t"
                f"{conc_pct_str}\t{total_gt_str}\t{concord_gt_str}\t{passes_str}\t"
                f"{r['BETA_pan']}\t{r['SE_pan']}\t{r['p_pan']}\t{r['AF_pan']}\t"
                f"{r['BETA_lin']}\t{r['SE_lin']}\t{r['p_lin']}\t{r['AF_lin']}\t"
                f"{db}\t{zb}\t{dlp}\t{zlp}\t{daf}\t{zaf}\n"
            )

    print(f"Wrote concordance-annotated pre-concordance metrics table: {out_tsv}")
    print(f"Post-concordance shared variant count (downstream set): {len(kept):,}")


if __name__ == "__main__":
    main()