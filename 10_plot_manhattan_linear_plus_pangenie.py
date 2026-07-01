#!/usr/bin/env python3
import argparse
import gzip
import math
import sys
from pathlib import Path
from collections import defaultdict
from statistics import median

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MIN_P = 1e-300

# Colorblind-friendly palette
COLOR_LINEAR = "0.55"
COLOR_UNIQUE = "#CC79A7"
COLOR_BOOSTED = "#009E73"


# -------------------------
# Helpers
# -------------------------
def open_maybe_gzip(path, mode="rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def ensure_parent_dir(pathlike):
    p = Path(pathlike)
    if p.parent and str(p.parent) not in (".", ""):
        p.parent.mkdir(parents=True, exist_ok=True)


def save_fig_pdf_png(fig, out_pdf, *, png_dpi=300):
    """
    Save one figure to both:
      - PDF at out_pdf
      - PNG at out_pdf with '.png' extension
    """
    ensure_parent_dir(out_pdf)
    out_pdf = str(out_pdf)
    if not out_pdf.lower().endswith(".pdf"):
        raise ValueError(f"Expected .pdf output path, got: {out_pdf}")
    out_png = out_pdf[:-4] + ".png"

    fig.savefig(out_pdf, format="pdf", bbox_inches="tight")
    fig.savefig(out_png, format="png", dpi=png_dpi, bbox_inches="tight")


def canonical_chrom(chrom: str) -> str:
    if chrom.lower().startswith("chr"):
        return chrom[3:]
    return chrom.lower()


def parse_marker_id(marker_id: str):
    parts = marker_id.split("_")
    if len(parts) < 2:
        raise ValueError
    chrom_raw = parts[0]
    pos1 = int(parts[1])
    ck = canonical_chrom(chrom_raw)
    return chrom_raw, pos1, ck


def safe_p(p):
    if p is None or not math.isfinite(p):
        return float("nan")
    if p <= 0:
        return MIN_P
    if p > 1:
        return 1.0
    return p


def safe_logp(p):
    p = safe_p(p)
    if not math.isfinite(p):
        return float("nan")
    return -math.log10(p)


def load_id_list(path):
    """
    Load MarkerIDs from:
      - plain list (one per line), OR
      - whitespace-delimited table containing a 'MarkerID' column
    """
    ids = set()
    with open_maybe_gzip(path, "rt") as f:
        first = f.readline().strip()
        if not first:
            return ids
        fields = first.split()
        if "MarkerID" in fields:
            mid_i = fields.index("MarkerID")
            for line in f:
                if line.strip():
                    parts = line.split()
                    if len(parts) > mid_i:
                        ids.add(parts[mid_i])
        else:
            ids.add(first)
            for line in f:
                if line.strip():
                    ids.add(line.strip())
    return ids


def parse_header_indices(header, path):
    cols = header.split()
    try:
        return cols.index("MarkerID"), cols.index("p.value")
    except ValueError:
        sys.exit(f"ERROR: MarkerID / p.value missing in {path}")


def chrom_sort_key(ck: str):
    """
    Sort chromosomes in natural order:
      1,2,...,22,X,Y,M, then others
    """
    c = ck.lower()
    try:
        return (0, int(c))
    except ValueError:
        pass
    if c == "x":
        return (1, 23)
    if c == "y":
        return (1, 24)
    if c in ("m", "mt"):
        return (1, 25)
    return (2, c)


def add_hline_label(ax, y, label, *, color="0.1"):
    """
    Add a text label at the left edge of the plot for a horizontal line.
    y is in data coords; x is in axes coords.
    """
    ax.text(
        0.01,
        y,
        label,
        transform=ax.get_yaxis_transform(),
        ha="left",
        va="bottom",
        fontsize=8,
        color=color,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.75, pad=1.5),
        zorder=10,
    )


def format_p(p: float) -> str:
    if p == 0:
        return "0"
    if p < 1e-3:
        return f"{p:.1e}"
    return f"{p:g}"


# -------------------------
# NEW: write boosted (pangenie-stronger/weaker) shared variants that are plotted
# -------------------------
def write_boosted_outputs(prefix, direction, boosted_records, nominal_p, p_fold, p_max):
    """
    boosted_records: list of dicts with keys:
      mid, ck, pos1, pP, logpP, pL, logpL, dlogp, fold

    Writes (exactly the boosted shared variants used for plotting, i.e. among pangenie p<=p_max):
      - {prefix}.boosted_shared.<direction>.ALL.ids
      - {prefix}.boosted_shared.<direction>.ALL.tsv
      - {prefix}.boosted_shared.<direction>.NOMINAL.ids
      - {prefix}.boosted_shared.<direction>.NOMINAL.tsv
    """
    def sort_key(r):
        return (chrom_sort_key(r["ck"]), r["pos1"], r["mid"])

    boosted_sorted = sorted(boosted_records, key=sort_key)
    boosted_nom = [r for r in boosted_sorted if r["pP"] <= nominal_p]

    tag = f"{prefix}.boosted_shared.{direction}"
    out_ids_all = f"{tag}.ALL.ids"
    out_tsv_all = f"{tag}.ALL.tsv"
    out_ids_nom = f"{tag}.NOMINAL.ids"
    out_tsv_nom = f"{tag}.NOMINAL.tsv"

    # IDs (ALL)
    ensure_parent_dir(out_ids_all)
    with open(out_ids_all, "w") as f:
        for r in boosted_sorted:
            f.write(r["mid"] + "\n")

    # TSV (ALL)
    ensure_parent_dir(out_tsv_all)
    with open(out_tsv_all, "w") as f:
        f.write(
            "MarkerID\tChrom\tPos\tpPangenie\tpLinear\tlogpPangenie\tlogpLinear\t"
            "dLogp(logpP-logpL)\tFold(pL/pP)\tIsNominalP\tNote\n"
        )
        note = f"boosted_shared among pangenie p<= {p_max:g}; criterion uses p_fold={p_fold:g} and direction={direction}"
        for r in boosted_sorted:
            f.write(
                f"{r['mid']}\t{r['ck']}\t{r['pos1']}\t{r['pP']:.12g}\t{r['pL']:.12g}\t"
                f"{r['logpP']:.6f}\t{r['logpL']:.6f}\t{r['dlogp']:.6f}\t{r['fold']:.6f}\t"
                f"{1 if r['pP'] <= nominal_p else 0}\t{note}\n"
            )

    # IDs (NOMINAL)
    ensure_parent_dir(out_ids_nom)
    with open(out_ids_nom, "w") as f:
        for r in boosted_nom:
            f.write(r["mid"] + "\n")

    # TSV (NOMINAL)
    ensure_parent_dir(out_tsv_nom)
    with open(out_tsv_nom, "w") as f:
        f.write(
            "MarkerID\tChrom\tPos\tpPangenie\tpLinear\tlogpPangenie\tlogpLinear\t"
            "dLogp(logpP-logpL)\tFold(pL/pP)\tIsNominalP\tNote\n"
        )
        note = f"boosted_shared (NOMINAL pP<= {nominal_p:g}) among pangenie p<= {p_max:g}; p_fold={p_fold:g}; direction={direction}"
        for r in boosted_nom:
            f.write(
                f"{r['mid']}\t{r['ck']}\t{r['pos1']}\t{r['pP']:.12g}\t{r['pL']:.12g}\t"
                f"{r['logpP']:.6f}\t{r['logpL']:.6f}\t{r['dlogp']:.6f}\t{r['fold']:.6f}\t1\t{note}\n"
            )

    print(f"Wrote boosted shared IDs (ALL):       {out_ids_all}  ({len(boosted_sorted):,})")
    print(f"Wrote boosted shared table (ALL):     {out_tsv_all}  ({len(boosted_sorted):,})")
    print(f"Wrote boosted shared IDs (NOMINAL):   {out_ids_nom}  ({len(boosted_nom):,})")
    print(f"Wrote boosted shared table (NOMINAL): {out_tsv_nom}  ({len(boosted_nom):,})")


# -------------------------
# QQ + Lambda helpers (no SciPy)
# -------------------------
def _norm_ppf(p: float) -> float:
    """
    Inverse CDF (ppf) for standard normal using Acklam's approximation.
    Accurate enough for QQ/lambda work; avoids SciPy dependency.
    """
    if p <= 0.0:
        return -math.inf
    if p >= 1.0:
        return math.inf

    # Coefficients (Peter John Acklam)
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00]

    plow = 0.02425
    phigh = 1 - plow

    if p < plow:
        q = math.sqrt(-2 * math.log(p))
        return (((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) / \
               ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1)
    if p > phigh:
        q = math.sqrt(-2 * math.log(1 - p))
        return -(((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) / \
                 ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1)

    q = p - 0.5
    r = q * q
    return (((((a[0]*r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5]) * q / \
           (((((b[0]*r + b[1])*r + b[2])*r + b[3])*r + b[4])*r + 1)


def chi2_1df_isf(p: float) -> float:
    """
    Inverse survival function for chi-square with 1 df.
    For 1 df, chi2 = z^2 where z = N(0,1) quantile for p/2 (two-sided).
    """
    p = safe_p(p)
    if not math.isfinite(p):
        return float("nan")
    if p >= 1.0:
        return 0.0
    if p <= 0.0:
        p = MIN_P

    p2 = p / 2.0
    if p2 < MIN_P:
        p2 = MIN_P
    if p2 > 0.5:
        p2 = 0.5

    z = _norm_ppf(1.0 - p2)
    return z * z


def lambda_gc_from_pvalues(pvals):
    """
    Genomic inflation factor lambda_GC:
      lambda = median(chi2_1df) / 0.4549364...
    """
    chi2_vals = []
    for p in pvals:
        x = chi2_1df_isf(p)
        if math.isfinite(x):
            chi2_vals.append(x)

    if not chi2_vals:
        return float("nan"), float("nan"), 0

    med = median(chi2_vals)
    chi2_median_null = 0.454936423119572
    lam = med / chi2_median_null
    return lam, med, len(chi2_vals)


def qq_plot(pvals, out_pdf, *, title, png_dpi=300, max_points=500_000):
    """
    Standard QQ plot: expected -log10(p) vs observed -log10(p).
    Downsamples plotting if pvals is huge (lambda should be computed separately on full set).
    """
    p_clean = [safe_p(p) for p in pvals if p is not None and math.isfinite(safe_p(p))]
    if not p_clean:
        print("WARNING: No valid p-values for QQ plot.")
        return

    p_sorted = sorted(p_clean)
    n = len(p_sorted)

    if n > max_points:
        step = max(1, n // max_points)
        idxs = list(range(step - 1, n, step))
        if idxs[-1] != n - 1:
            idxs.append(n - 1)
    else:
        idxs = list(range(n))

    exp = []
    obs = []
    for i in idxs:
        rank = i + 1
        p_exp = rank / (n + 1.0)
        exp.append(-math.log10(p_exp))
        obs.append(-math.log10(max(p_sorted[i], MIN_P)))

    fig, ax = plt.subplots(figsize=(6, 6), dpi=200)
    ax.scatter(exp, obs, s=6, marker=".", linewidths=0, alpha=0.8, rasterized=True)

    mx = max(exp[-1], obs[-1]) if exp and obs else 1.0
    

    ax.set_xlabel("Expected -log10(p)")
    ax.set_ylabel("Observed -log10(p)")
    ax.set_title(title)

    fig.tight_layout()
    save_fig_pdf_png(fig, out_pdf, png_dpi=png_dpi)
    plt.close(fig)
    print(f"Wrote: {out_pdf} and {out_pdf[:-4] + '.png'}")


# -------------------------
# NEW: nominal TSV writer
# -------------------------
def pad_fields(fields, ncols):
    if len(fields) < ncols:
        return fields + [""] * (ncols - len(fields))
    return fields[:ncols]


def write_labeled_nominal_tsv(
    out_tsv,
    *,
    nominal_p,
    hap_unique,
    boosted_ids,
    pangenie_header,
    pangenie_rows,   # mid -> fields list (full row)
    linear_header,
    linear_rows,     # mid -> fields list (full row)
):
    """
    Writes one TSV for all nominal hits, labeled as boosted/unique/linear.
    Includes both GWAS rows (prefixed P_ and L_) to preserve "GWAS statistics".
    """
    ensure_parent_dir(out_tsv)

    # Build union of nominal IDs captured from either file
    mids = set(pangenie_rows.keys()) | set(linear_rows.keys())
    if not mids:
        print("No nominal variants captured; skipping labeled nominal TSV.")
        return

    # Identify p.value indices (needed for pP/pL columns)
    def idx(cols, name):
        try:
            return cols.index(name)
        except ValueError:
            return None

    p_p_i = idx(pangenie_header, "p.value")
    l_p_i = idx(linear_header, "p.value")

    # Sort by chrom/pos when possible
    sortable = []
    unsortable = []
    for mid in mids:
        try:
            _, pos1, ck = parse_marker_id(mid)
            sortable.append((chrom_sort_key(ck), pos1, mid, ck, pos1))
        except Exception:
            unsortable.append(mid)

    sortable.sort()
    ordered = [mid for _, _, mid, _, _ in sortable] + sorted(unsortable)

    # Output header
    # Keep a compact computed-stat block + all original columns (prefixed)
    base_cols = [
        "Label",
        "MarkerID",
        "Chrom",
        "Pos",
        "IsNominalPangenie",
        "IsNominalLinear",
        "pPangenie",
        "pLinear",
        "logpPangenie",
        "logpLinear",
        "dLogp(logpP-logpL)",
        "Fold(pL/pP)",
    ]

    # Prefix original GWAS columns to avoid collisions
    p_pref_cols = ["P_" + c for c in pangenie_header if c != "MarkerID"]
    l_pref_cols = ["L_" + c for c in linear_header if c != "MarkerID"]

    with open(out_tsv, "w") as out:
        out.write("\t".join(base_cols + p_pref_cols + l_pref_cols) + "\n")

        for mid in ordered:
            # chrom/pos best-effort
            ck = ""
            pos1 = ""
            try:
                _, pos_int, ck_int = parse_marker_id(mid)
                ck = ck_int
                pos1 = str(pos_int)
            except Exception:
                pass

            # label priority: boosted > unique > linear
            if mid in boosted_ids:
                label = "boosted"
            elif mid in hap_unique:
                label = "unique"
            else:
                label = "linear"

            # pull p-values if present
            pP = float("nan")
            pL = float("nan")
            lpP = float("nan")
            lpL = float("nan")

            p_fields = pangenie_rows.get(mid)
            l_fields = linear_rows.get(mid)

            is_nom_p = 0
            is_nom_l = 0

            if p_fields is not None and p_p_i is not None and p_p_i < len(p_fields):
                try:
                    pP = safe_p(float(p_fields[p_p_i]))
                except Exception:
                    pP = float("nan")
                if math.isfinite(pP):
                    lpP = safe_logp(pP)
                    if pP <= nominal_p:
                        is_nom_p = 1

            if l_fields is not None and l_p_i is not None and l_p_i < len(l_fields):
                try:
                    pL = safe_p(float(l_fields[l_p_i]))
                except Exception:
                    pL = float("nan")
                if math.isfinite(pL):
                    lpL = safe_logp(pL)
                    if pL <= nominal_p:
                        is_nom_l = 1

            # derived
            dlogp = float("nan")
            fold = float("nan")
            if math.isfinite(lpP) and math.isfinite(lpL):
                dlogp = lpP - lpL
            if math.isfinite(pP) and math.isfinite(pL) and pP > 0:
                fold = pL / pP

            def fmt(x, nd=6):
                if x is None or not math.isfinite(x):
                    return "NA"
                return f"{x:.{nd}f}"

            base_vals = [
                label,
                mid,
                ck if ck else "NA",
                pos1 if pos1 else "NA",
                str(is_nom_p),
                str(is_nom_l),
                (f"{pP:.12g}" if math.isfinite(pP) else "NA"),
                (f"{pL:.12g}" if math.isfinite(pL) else "NA"),
                fmt(lpP, 6),
                fmt(lpL, 6),
                fmt(dlogp, 6),
                fmt(fold, 6),
            ]

            # prefixed row expansions (skip MarkerID)
            p_out = []
            if p_fields is not None:
                p_fields = pad_fields(p_fields, len(pangenie_header))
                for j, c in enumerate(pangenie_header):
                    if c == "MarkerID":
                        continue
                    p_out.append(p_fields[j])
            else:
                p_out = [""] * len(p_pref_cols)

            l_out = []
            if l_fields is not None:
                l_fields = pad_fields(l_fields, len(linear_header))
                for j, c in enumerate(linear_header):
                    if c == "MarkerID":
                        continue
                    l_out.append(l_fields[j])
            else:
                l_out = [""] * len(l_pref_cols)

            out.write("\t".join(base_vals + p_out + l_out) + "\n")

    print(f"Wrote labeled nominal variants TSV: {out_tsv}  ({len(mids):,} variants)")


# -------------------------
# Main
# -------------------------
def main():
    matplotlib.rcParams["agg.path.chunksize"] = 20000

    ap = argparse.ArgumentParser(
        description="Manhattan plot using hap.py–defined unique and shared variants (PDF+PNG outputs)"
    )
    ap.add_argument("linear_gwas")
    ap.add_argument("pangenie_gwas")
    ap.add_argument("--shared-ids", required=True)
    ap.add_argument("--pangenie-unique-ids", required=True)
    ap.add_argument("--prefix", default="manhattan")

    ap.add_argument("--p-max", type=float, default=0.01)
    ap.add_argument("--nominal-p", type=float, default=5e-5)
    ap.add_argument("--gwas-p", type=float, default=5e-8)
    ap.add_argument("--p-fold", type=float, default=10.0)
    ap.add_argument(
        "--direction",
        choices=["pangenie-stronger", "pangenie-weaker", "both"],
        default="pangenie-stronger",
    )

    ap.add_argument("--linear-point-size", type=float, default=6.0)
    ap.add_argument("--unique-point-size", type=float, default=10.0)
    ap.add_argument("--boosted-point-size", type=float, default=10.0)
    ap.add_argument("--alpha-linear", type=float, default=1.0)
    ap.add_argument("--alpha-unique", type=float, default=0.9)
    ap.add_argument("--alpha-boosted", type=float, default=0.95)

    ap.add_argument("--dpi", type=int, default=200)
    ap.add_argument("--png-dpi", type=int, default=300, help="PNG DPI (default 300)")

    ap.add_argument(
        "--exclude-ids",
        default=None,
        help="Optional file containing MarkerIDs to exclude (one per line or table w/ MarkerID column).",
    )
    ap.add_argument(
        "--exclude-id",
        action="append",
        default=[],
        help="MarkerID to exclude (repeatable). Example: --exclude-id chr1_123_A_G",
    )

    # QQ options
    ap.add_argument("--no-qq", action="store_true", help="Disable Pangenie QQ plot + lambda")
    ap.add_argument("--qq-max-points", type=int, default=500_000, help="Max points to plot in QQ (downsample if larger)")

    # boosted outputs toggle
    ap.add_argument(
        "--no-write-boosted",
        action="store_true",
        help="Disable writing boosted shared outputs (ids/tsv).",
    )

    # NEW: nominal TSV output
    ap.add_argument(
        "--no-nominal-tsv",
        action="store_true",
        help="Disable writing a labeled nominal variants TSV.",
    )
    ap.add_argument(
        "--nominal-tsv",
        default=None,
        help="Output path for labeled nominal variants TSV (default: <prefix>.nominal_variants.labeled.tsv).",
    )

    args = ap.parse_args()

    delta_thr = math.log10(args.p_fold)

    # hap.py-derived sets
    hap_shared = load_id_list(args.shared_ids)
    hap_unique = load_id_list(args.pangenie_unique_ids)

    # Exclusions
    exclude = set()
    if args.exclude_ids:
        exclude |= load_id_list(args.exclude_ids)
    exclude |= set(args.exclude_id)

    if exclude:
        hap_shared -= exclude
        hap_unique -= exclude

    print(f"hap.py shared IDs:   {len(hap_shared):,}")
    print(f"hap.py unique IDs:   {len(hap_unique):,}")
    print(f"excluded IDs:        {len(exclude):,}")

    # -------------------------
    # Load Pangenie candidates (p <= p-max)
    # ALSO collect all p-values for QQ/lambda (unfiltered by p-max)
    # ALSO collect nominal rows for TSV export
    # -------------------------
    p_cand = {}  # mid -> (ck, pos1, pP, logpP)
    chrom_max = defaultdict(int)
    pangenie_pvals_all = []

    pangenie_header = None
    pangenie_nominal_rows = {}  # mid -> full fields list

    with open_maybe_gzip(args.pangenie_gwas, "rt") as f:
        header = f.readline().strip()
        pangenie_header = header.split()
        mid_i, p_i = parse_header_indices(header, args.pangenie_gwas)

        for line in f:
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) <= max(mid_i, p_i):
                continue

            mid = fields[mid_i]
            if mid in exclude:
                continue

            try:
                pP = float(fields[p_i])
            except ValueError:
                continue
            pP = safe_p(pP)
            if not math.isfinite(pP):
                continue

            pangenie_pvals_all.append(pP)

            # store nominal rows (for output)
            if pP <= args.nominal_p:
                pangenie_nominal_rows[mid] = pad_fields(fields, len(pangenie_header))

            try:
                _, pos1, ck = parse_marker_id(mid)
            except Exception:
                continue

            chrom_max[ck] = max(chrom_max[ck], pos1)

            if pP <= args.p_max:
                p_cand[mid] = (ck, pos1, pP, safe_logp(pP))

    # -------------------------
    # Load linear p-values for hap.py shared variants (for boost calc)
    # -------------------------
    linear_stats = {}  # mid -> (pL, logpL)

    with open_maybe_gzip(args.linear_gwas, "rt") as f:
        header = f.readline().strip()
        mid_i, p_i = parse_header_indices(header, args.linear_gwas)

        for line in f:
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) <= max(mid_i, p_i):
                continue

            mid = fields[mid_i]
            if mid in exclude:
                continue
            if mid not in hap_shared:
                continue

            try:
                pL = float(fields[p_i])
            except ValueError:
                continue
            pL = safe_p(pL)
            if math.isfinite(pL):
                linear_stats[mid] = (pL, safe_logp(pL))

    # -------------------------
    # Build plot lists for overlays (+ boosted output records)
    # -------------------------
    unique_all = []
    unique_nom = []
    boosted_all = []
    boosted_nom = []
    boosted_records = []

    for mid, (ck, pos1, pP, logpP) in p_cand.items():
        if mid in hap_unique:
            unique_all.append((ck, pos1, logpP))
            if pP <= args.nominal_p:
                unique_nom.append((ck, pos1, logpP))
        elif mid in hap_shared and mid in linear_stats:
            pL, logpL = linear_stats[mid]
            dlogp = logpP - logpL  # = log10(pL/pP)

            keep = (
                (args.direction == "pangenie-stronger" and dlogp >= delta_thr) or
                (args.direction == "pangenie-weaker" and dlogp <= -delta_thr) or
                (args.direction == "both" and abs(dlogp) >= delta_thr)
            )
            if keep:
                boosted_all.append((ck, pos1, logpP))
                boosted_records.append(
                    {
                        "mid": mid,
                        "ck": ck,
                        "pos1": pos1,
                        "pP": pP,
                        "logpP": logpP,
                        "pL": pL,
                        "logpL": logpL,
                        "dlogp": dlogp,
                        "fold": 10 ** dlogp,  # pL/pP
                    }
                )
                if pP <= args.nominal_p:
                    boosted_nom.append((ck, pos1, logpP))

    boosted_ids = set(r["mid"] for r in boosted_records)

    print(f"Pangenie unique (p<=p-max): {len(unique_all):,}")
    print(f"Boosted shared (p-fold):    {len(boosted_all):,}")

    if not args.no_write_boosted:
        write_boosted_outputs(
            args.prefix,
            args.direction,
            boosted_records,
            nominal_p=args.nominal_p,
            p_fold=args.p_fold,
            p_max=args.p_max,
        )

    # =========================
    # Pangenie QQ plot + lambda_GC
    # =========================
    if not args.no_qq:
        lam, med_chi2, n_used = lambda_gc_from_pvalues(pangenie_pvals_all)
        if math.isfinite(lam):
            print(f"Pangenie lambda_GC: {lam:.4f} (n={n_used:,}, median_chi2={med_chi2:.4f})")
        else:
            print("Pangenie lambda_GC: NA (no valid p-values)")

        qq_out = f"{args.prefix}.qq.pangenie.pdf"
        title = f"Pangenie GWAS QQ (λGC={lam:.3f})" if math.isfinite(lam) else "Pangenie GWAS QQ"
        qq_plot(
            pangenie_pvals_all,
            qq_out,
            title=title,
            png_dpi=args.png_dpi,
            max_points=args.qq_max_points,
        )

    # =========================
    # Manhattan plotting section
    # =========================
    if not chrom_max:
        sys.exit("ERROR: No chromosome positions found from Pangenie GWAS input after filtering.")

    chroms = sorted(chrom_max.keys(), key=chrom_sort_key)
    offsets = {}
    centers = {}
    cum = 0
    gap = 5_000_000
    for ck in chroms:
        offsets[ck] = cum
        length = chrom_max[ck]
        centers[ck] = cum + length / 2
        cum += length + gap

    def to_x(ck, pos1):
        return offsets[ck] + pos1

    y_gwas = safe_logp(args.gwas_p)
    y_nom = safe_logp(args.nominal_p)

    all_y = []
    all_y.extend([lp for _, _, lp in unique_all])
    all_y.extend([lp for _, _, lp in boosted_all])
    max_y = max(all_y + [y_gwas, y_nom, 10.0])
    y_top = max_y * 1.1

    fig, ax = plt.subplots(figsize=(14, 6), dpi=args.dpi)

    # ---- Linear background ----
    xs, ys = [], []

    # NEW: capture nominal linear rows for TSV export while we stream the file
    linear_header = None
    linear_nominal_rows = {}

    with open_maybe_gzip(args.linear_gwas, "rt") as f:
        header = f.readline().strip()
        linear_header = header.split()
        mid_i, p_i = parse_header_indices(header, args.linear_gwas)

        for line in f:
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) <= max(mid_i, p_i):
                continue

            mid = fields[mid_i]
            if mid in exclude:
                continue

            try:
                pL = float(fields[p_i])
            except ValueError:
                continue
            pL = safe_p(pL)
            if not math.isfinite(pL) or pL > args.p_max:
                continue

            # store nominal rows (for output)
            if pL <= args.nominal_p:
                linear_nominal_rows[mid] = pad_fields(fields, len(linear_header))

            try:
                _, pos1, ck = parse_marker_id(mid)
            except Exception:
                continue

            if ck not in offsets:
                continue

            xs.append(to_x(ck, pos1))
            ys.append(safe_logp(pL))

    ax.scatter(
        xs,
        ys,
        s=args.linear_point_size,
        alpha=args.alpha_linear,
        color=COLOR_LINEAR,
        linewidths=0,
        marker=".",
        zorder=2,
        rasterized=True,
        label=f"Linear (p≤{args.p_max:g})",
    )

    ax.axhline(y_gwas, linestyle="--", color="0.2", linewidth=1)
    ax.axhline(y_nom, linestyle="--", color="0.2", linewidth=1)
    add_hline_label(ax, y_gwas, f"GWAS: p={format_p(args.gwas_p)}")
    add_hline_label(ax, y_nom, f"Nominal: p={format_p(args.nominal_p)}")

    ax.set_xticks([centers[ck] for ck in chroms])
    ax.set_xticklabels([ck.upper() if not ck.isdigit() else ck for ck in chroms], fontsize=8)

    ax.set_ylim(0, y_top)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("-log10(p.value)")

    ax.scatter([], [], color=COLOR_LINEAR, s=args.linear_point_size, label=f"Linear (p≤{args.p_max:g})")

    # NEW: write labeled nominal TSV (after we've collected both pangenie + linear nominal rows)
    if not args.no_nominal_tsv:
        out_tsv = args.nominal_tsv or f"{args.prefix}.nominal_variants.labeled.tsv"
        write_labeled_nominal_tsv(
            out_tsv,
            nominal_p=args.nominal_p,
            hap_unique=hap_unique,
            boosted_ids=boosted_ids,
            pangenie_header=pangenie_header,
            pangenie_rows=pangenie_nominal_rows,
            linear_header=linear_header,
            linear_rows=linear_nominal_rows,
        )

    def overlay_and_save(unique_pts, boosted_pts, title, out_pdf):
        ensure_parent_dir(out_pdf)
        artists = []

        if unique_pts:
            xu = [to_x(ck, pos1) for ck, pos1, _ in unique_pts]
            yu = [lp for _, _, lp in unique_pts]
            sc_u = ax.scatter(
                xu,
                yu,
                s=args.unique_point_size,
                alpha=args.alpha_unique,
                color=COLOR_UNIQUE,
                linewidths=0,
                zorder=4,
                label="Pangenie-unique",
            )
            artists.append(sc_u)

        if boosted_pts:
            xb = [to_x(ck, pos1) for ck, pos1, _ in boosted_pts]
            yb = [lp for _, _, lp in boosted_pts]
            sc_b = ax.scatter(
                xb,
                yb,
                s=args.boosted_point_size,
                alpha=args.alpha_boosted,
                color=COLOR_BOOSTED,
                linewidths=0,
                zorder=5,
                label=f"Boosted shared (fold≥{args.p_fold:g}, {args.direction})",
            )
            artists.append(sc_b)

        ax.set_title(title)
        leg = ax.legend(loc="upper right", fontsize=8, frameon=True)

        fig.tight_layout()
        save_fig_pdf_png(fig, out_pdf, png_dpi=args.png_dpi)
        print(f"Wrote: {out_pdf} and {out_pdf[:-4] + '.png'}")

        leg.remove()
        for a in artists:
            a.remove()

    overlay_and_save(
        unique_all,
        boosted_all,
        title="Manhattan: Linear background + ALL Pangenie-unique + ALL boosted shared",
        out_pdf=f"{args.prefix}.manhattan.unique_plus_boosted_ALL.pdf",
    )

    overlay_and_save(
        unique_nom,
        boosted_nom,
        title=f"Manhattan: Linear background + NOMINAL unique + NOMINAL boosted (p≤{format_p(args.nominal_p)})",
        out_pdf=f"{args.prefix}.manhattan.unique_plus_boosted_NOMINAL.pdf",
    )

    plt.close(fig)


if __name__ == "__main__":
    main()

