#!/usr/bin/env python3
import argparse
import gzip
import sys
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Set, Tuple, Optional


def open_maybe_gz(path: str, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def norm(s: str) -> str:
    return s.strip().strip('"').replace("\r", "")


def load_id_list(path: str) -> Set[str]:
    """
    Accepts:
      - one ID per line, OR
      - whitespace/tab delimited table with header containing 'MarkerID'
    Returns set of MarkerIDs.
    """
    ids: Set[str] = set()
    with open_maybe_gz(path, "rt") as f:
        # find first non-empty line
        first = None
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            first = line
            break
        if first is None:
            return ids

        fields = first.split()
        if "MarkerID" in fields:
            mid_i = fields.index("MarkerID")
            for line in f:
                if not line.strip() or line.lstrip().startswith("#"):
                    continue
                parts = line.split()
                if len(parts) > mid_i:
                    ids.add(norm(parts[mid_i]))
        else:
            # first line is an ID
            ids.add(norm(fields[0]))
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                ids.add(norm(line.split()[0]))
    return ids


def load_boosted_ids(path: str) -> Set[str]:
    """One ID per line (first field), ignore comments/blank."""
    out = set()
    with open_maybe_gz(path, "rt") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            out.add(norm(line.split()[0]))
    return out


def parse_leads_table(leads_tsv: str) -> Tuple[List[str], Dict[str, Dict[str, str]]]:
    """
    leads.from_clumped.tsv has at least columns: CHR SNP BP P TOTAL NSIG (names vary slightly)
    Returns:
      - lead order list
      - lead_info dict: lead -> {CHR, BP, P}
    """
    with open_maybe_gz(leads_tsv, "rt") as f:
        header = f.readline().strip().split()
        if not header:
            raise RuntimeError(f"Empty leads table: {leads_tsv}")

        def find_col(*names):
            for n in names:
                if n in header:
                    return header.index(n)
            return None

        i_chr = find_col("CHR") or 0
        i_snp = find_col("SNP") or 1
        i_bp = find_col("BP") or 2
        i_p = find_col("P_CLUMP", "P", "LEAD_P") or 3

        order: List[str] = []
        info: Dict[str, Dict[str, str]] = {}

        for line in f:
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) <= max(i_chr, i_snp, i_bp, i_p):
                continue
            lead = norm(parts[i_snp])
            if not lead:
                continue
            order.append(lead)
            info[lead] = {
                "CHR": parts[i_chr],
                "BP": parts[i_bp],
                "P": parts[i_p],
            }

    return order, info


def parse_clumped_members(clumped_path: str) -> Tuple[Dict[str, List[str]], Set[str]]:
    """
    Parse PLINK .clumped:
      - lead in SNP column
      - members in SP2 column (comma-separated, with '(p)' suffix)
    Returns:
      members_by_lead: lead -> list of member IDs (including lead)
      all_members: set of all member IDs seen
    """
    members_by_lead: Dict[str, List[str]] = {}
    all_members: Set[str] = set()

    with open_maybe_gz(clumped_path, "rt") as f:
        header_line = f.readline().strip()
        if not header_line:
            raise RuntimeError(f"Empty clumped file: {clumped_path}")
        header = header_line.split()
        # typical columns: CHR F SNP BP P TOTAL NSIG ... SP2
        i_snp = header.index("SNP") if "SNP" in header else 2
        i_sp2 = header.index("SP2") if "SP2" in header else None

        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) <= i_snp:
                continue
            lead = norm(parts[i_snp])
            if not lead:
                continue

            mems = [lead]
            if i_sp2 is not None and i_sp2 < len(parts):
                sp2 = norm(parts[i_sp2])
                if sp2 and sp2 not in ("NONE", "NA"):
                    for tok in sp2.split(","):
                        tok = tok.strip()
                        if not tok:
                            continue
                        # strip "(p)"
                        if "(" in tok:
                            tok = tok.split("(", 1)[0].strip()
                        tok = norm(tok)
                        if tok and tok != "NONE":
                            mems.append(tok)

            # dedupe while preserving order
            seen = set()
            deduped = []
            for m in mems:
                if m not in seen:
                    seen.add(m)
                    deduped.append(m)

            members_by_lead[lead] = deduped
            all_members.update(deduped)

    return members_by_lead, all_members


def load_p_clump_for_members(clump_table: str, member_ids: Set[str]) -> Dict[str, float]:
    """
    clump table: OUT.for_clump.tsv with header containing MarkerID and P_CLUMP
    Only store p-values for member_ids to save memory.
    """
    pvals: Dict[str, float] = {}
    with open_maybe_gz(clump_table, "rt") as f:
        header = f.readline().strip().split()
        if not header:
            raise RuntimeError(f"Empty clump table: {clump_table}")
        if "MarkerID" not in header or "P_CLUMP" not in header:
            raise RuntimeError("clump-table must contain MarkerID and P_CLUMP columns")
        i_mid = header.index("MarkerID")
        i_p = header.index("P_CLUMP")

        for line in f:
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) <= max(i_mid, i_p):
                continue
            mid = norm(parts[i_mid])
            if mid not in member_ids:
                continue
            try:
                p = float(parts[i_p])
            except ValueError:
                continue
            if p <= 0:
                p = 1e-300
            pvals[mid] = p
    return pvals


def compute_nominal_member_counts(
    members_by_lead: Dict[str, List[str]],
    p_clump: Dict[str, float],
    override: Set[str],
    p_thr: float,
) -> Dict[str, Tuple[int, int, int]]:
    """
    Returns lead -> (nom_total, nom_override, nom_nonoverride)
    """
    out: Dict[str, Tuple[int, int, int]] = {}
    for lead, mems in members_by_lead.items():
        nt = no = nn = 0
        for m in mems:
            p = p_clump.get(m)
            if p is None:
                continue
            if p < p_thr:
                nt += 1
                if m in override:
                    no += 1
                else:
                    nn += 1
        out[lead] = (nt, no, nn)
    return out


def parse_ld_counts_and_max(
    ld_file: str,
    leads: Set[str],
    unique: Set[str],
    boosted: Set[str],
    r2_link: float,
) -> Tuple[Dict[str, Set[str]], Dict[str, Set[str]], Dict[str, float], Dict[str, float]]:
    """
    For each lead SNP_A:
      - unique_hits[lead] = set of SNP_B in unique with R2>=r2_link (exclude self)
      - boosted_hits[lead] = set of SNP_B in boosted with R2>=r2_link (exclude self)
      - maxU[lead] = max R2 to any unique SNP_B (exclude self)
      - maxB[lead] = max R2 to any boosted SNP_B (exclude self)
    """
    unique_hits: Dict[str, Set[str]] = defaultdict(set)
    boosted_hits: Dict[str, Set[str]] = defaultdict(set)
    maxU: Dict[str, float] = {}
    maxB: Dict[str, float] = {}

    with open_maybe_gz(ld_file, "rt") as f:
        header = f.readline().strip().split()
        if not header:
            raise RuntimeError(f"Empty LD file: {ld_file}")

        # column positions by name if present
        def idx(name: str) -> Optional[int]:
            return header.index(name) if name in header else None

        i_a = idx("SNP_A")
        i_b = idx("SNP_B")
        i_r = idx("R2")
        if i_a is None or i_b is None or i_r is None:
            raise RuntimeError("LD header must contain SNP_A SNP_B R2")

        for line in f:
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) <= max(i_a, i_b, i_r):
                continue
            a = norm(parts[i_a])
            if a not in leads:
                continue
            b = norm(parts[i_b])
            if a == b:
                continue
            try:
                r = float(parts[i_r])
            except ValueError:
                continue

            if b in unique:
                prev = maxU.get(a, -1.0)
                if r > prev:
                    maxU[a] = r
                if r >= r2_link:
                    unique_hits[a].add(b)

            if b in boosted:
                prev = maxB.get(a, -1.0)
                if r > prev:
                    maxB[a] = r
                if r >= r2_link:
                    boosted_hits[a].add(b)

    return unique_hits, boosted_hits, maxU, maxB


def write_locus_stats(
    out_path: str,
    lead_order: List[str],
    lead_info: Dict[str, Dict[str, str]],
    unique: Set[str],
    boosted: Set[str],
    unique_hits: Dict[str, Set[str]],
    boosted_hits: Dict[str, Set[str]],
    nominal_counts: Dict[str, Tuple[int, int, int]],
):
    with open(out_path, "w") as out:
        out.write(
            "\t".join([
                "CHR", "BP", "LEAD_SNP", "LEAD_P",
                "LEAD_IS_UNIQUE", "LEAD_IS_BOOSTED",
                "UNIQUE_LD_GE_THR_EXCL_SELF", "BOOSTED_LD_GE_THR_EXCL_SELF",
                "UNIQUE_LD_GE_THR_INCL_SELF", "BOOSTED_LD_GE_THR_INCL_SELF",
                "HAS_UNIQUE_LD_GE_THR_INCL_SELF", "HAS_BOOSTED_LD_GE_THR_INCL_SELF", "HAS_ANY_LD_GE_THR_INCL_SELF",
                "NOMINAL_MEMBER_TOTAL", "NOMINAL_OVERRIDE_COUNT", "NOMINAL_NONOVERRIDE_COUNT",
                "IS_NOVEL_ONLY_OVERRIDE_NOMINAL", "NOVEL_CLASS_BY_LEAD",
            ]) + "\n"
        )

        for lead in lead_order:
            info = lead_info.get(lead, {"CHR": "NA", "BP": "NA", "P": "NA"})
            lead_is_u = 1 if lead in unique else 0
            lead_is_b = 1 if lead in boosted else 0

            u_ex = len(unique_hits.get(lead, set()))
            b_ex = len(boosted_hits.get(lead, set()))
            u_in = u_ex + lead_is_u
            b_in = b_ex + lead_is_b

            has_u = 1 if u_in > 0 else 0
            has_b = 1 if b_in > 0 else 0
            has_any = 1 if (has_u or has_b) else 0

            nt, no, nn = nominal_counts.get(lead, (0, 0, 0))
            is_novel = 1 if (nn == 0 and no > 0) else 0

            novel_class = "NA"
            if is_novel:
                if lead_is_u:
                    novel_class = "unique_lead"
                elif lead_is_b:
                    novel_class = "boosted_lead"
                else:
                    novel_class = "other_lead"

            out.write(
                "\t".join(map(str, [
                    info["CHR"], info["BP"], lead, info["P"],
                    lead_is_u, lead_is_b,
                    u_ex, b_ex,
                    u_in, b_in,
                    has_u, has_b, has_any,
                    nt, no, nn,
                    is_novel, novel_class,
                ])) + "\n"
            )


def write_overall_stats(
    out_path: str,
    locus_stats_path: str,
    r2_link: float,
    p_thr: float,
):
    # Aggregate from locus_stats table (simple and robust)
    total_loci = 0
    lead_unique = 0
    lead_boosted = 0
    loci_has_u = 0
    loci_has_b = 0
    loci_has_any = 0
    loci_has_both = 0
    sum_unique_in_ld = 0
    sum_boosted_in_ld = 0
    novel_total = 0
    novel_unique_lead = 0
    novel_boosted_lead = 0
    novel_other_lead = 0

    with open(locus_stats_path, "r") as f:
        header = f.readline().strip().split("\t")
        col = {name: i for i, name in enumerate(header)}

        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            total_loci += 1

            lead_is_u = int(parts[col["LEAD_IS_UNIQUE"]])
            lead_is_b = int(parts[col["LEAD_IS_BOOSTED"]])
            u_in = int(parts[col["UNIQUE_LD_GE_THR_INCL_SELF"]])
            b_in = int(parts[col["BOOSTED_LD_GE_THR_INCL_SELF"]])
            has_u = int(parts[col["HAS_UNIQUE_LD_GE_THR_INCL_SELF"]])
            has_b = int(parts[col["HAS_BOOSTED_LD_GE_THR_INCL_SELF"]])
            has_any = int(parts[col["HAS_ANY_LD_GE_THR_INCL_SELF"]])
            is_novel = int(parts[col["IS_NOVEL_ONLY_OVERRIDE_NOMINAL"]])
            nclass = parts[col["NOVEL_CLASS_BY_LEAD"]]

            if lead_is_u:
                lead_unique += 1
            if lead_is_b:
                lead_boosted += 1

            if has_u:
                loci_has_u += 1
                sum_unique_in_ld += u_in
            if has_b:
                loci_has_b += 1
                sum_boosted_in_ld += b_in
            if has_any:
                loci_has_any += 1
            if has_u and has_b:
                loci_has_both += 1

            if is_novel:
                novel_total += 1
                if nclass == "unique_lead":
                    novel_unique_lead += 1
                elif nclass == "boosted_lead":
                    novel_boosted_lead += 1
                else:
                    novel_other_lead += 1

    with open(out_path, "w") as out:
        def w(k, v): out.write(f"{k}\t{v}\n")

        w("TOTAL_NOMINAL_LOCI", total_loci)
        w("R2_LINK_THRESHOLD", r2_link)
        w("P_THRESHOLD_NOMINAL", p_thr)

        w("LEAD_IS_UNIQUE_COUNT", lead_unique)
        w("LEAD_IS_BOOSTED_COUNT", lead_boosted)

        w("LOCI_WITH_UNIQUE_IN_LD_INCL_SELF", loci_has_u)
        w("LOCI_WITH_BOOSTED_IN_LD_INCL_SELF", loci_has_b)
        w("LOCI_WITH_BOTH_UNIQUE_AND_BOOSTED_IN_LD_INCL_SELF", loci_has_both)
        w("LOCI_WITH_EITHER_UNIQUE_OR_BOOSTED_IN_LD_INCL_SELF", loci_has_any)

        w("SUM_DISTINCT_UNIQUE_VARIANTS_IN_LD_ACROSS_LOCI_INCL_SELF", sum_unique_in_ld)
        w("SUM_DISTINCT_BOOSTED_VARIANTS_IN_LD_ACROSS_LOCI_INCL_SELF", sum_boosted_in_ld)

        w("NOVEL_LOCI_ONLY_OVERRIDE_NOMINAL_TOTAL", novel_total)
        w("NOVEL_LOCI_ONLY_OVERRIDE_NOMINAL_LEAD_UNIQUE", novel_unique_lead)
        w("NOVEL_LOCI_ONLY_OVERRIDE_NOMINAL_LEAD_BOOSTED", novel_boosted_lead)
        if novel_other_lead > 0:
            w("NOVEL_LOCI_ONLY_OVERRIDE_NOMINAL_LEAD_OTHER", novel_other_lead)


def write_ld_sweep(
    out_path: str,
    lead_order: List[str],
    unique: Set[str],
    boosted: Set[str],
    maxU: Dict[str, float],
    maxB: Dict[str, float],
    sweep_min: float,
    sweep_max: float,
    sweep_step: float,
):
    # Build thresholds
    thr = []
    t = sweep_min
    # numeric stability
    while t <= sweep_max + 1e-12:
        thr.append(round(t, 2))
        t += sweep_step

    with open(out_path, "w") as out:
        out.write("\t".join([
            "R2_THRESHOLD",
            "N_LOCI_TOTAL",
            "N_LOCI_HAS_UNIQUE_LD_INCL_SELF",
            "N_LOCI_HAS_BOOSTED_LD_INCL_SELF",
            "N_LOCI_HAS_EITHER_INCL_SELF",
            "N_LOCI_HAS_BOTH_INCL_SELF",
        ]) + "\n")

        n_total = len(lead_order)
        for t in thr:
            hasU = hasB = hasEither = hasBoth = 0
            for lead in lead_order:
                u = (lead in unique) or (maxU.get(lead, -1.0) >= t)
                b = (lead in boosted) or (maxB.get(lead, -1.0) >= t)
                if u:
                    hasU += 1
                if b:
                    hasB += 1
                if u or b:
                    hasEither += 1
                if u and b:
                    hasBoth += 1
            out.write(f"{t:.2f}\t{n_total}\t{hasU}\t{hasB}\t{hasEither}\t{hasBoth}\n")

def classify_variant(marker_id: str) -> str:
    """Simple type classifier consistent with earlier scripts."""
    if len(marker_id) >= 50:
        return "SV"
    parts = marker_id.split("_")
    if len(parts) < 4:
        return "SNV"
    ref = parts[-2]
    alt = parts[-1]
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    return "INDEL"

def compute_variant_type_breakdown_from_ids(
    out_path: str,
    lead_order: List[str],
    unique_hits: Dict[str, Set[str]],
    boosted_hits: Dict[str, Set[str]],
    unique: Set[str],
    boosted: Set[str],
):
    """
    Combine unique + boosted variants in LD >= threshold (incl self if applicable),
    deduplicate across loci, then compute variant type breakdown
    using classify_variant().
    """

    combined_variants: Set[str] = set()

    for lead in lead_order:
        # LD hits (already filtered by r2_link upstream)
        combined_variants.update(unique_hits.get(lead, set()))
        combined_variants.update(boosted_hits.get(lead, set()))

        # include lead itself if it is unique or boosted
        if lead in unique or lead in boosted:
            combined_variants.add(lead)

    # Count types
    type_counts: Dict[str, int] = defaultdict(int)

    for vid in combined_variants:
        vtype = classify_variant(vid)
        type_counts[vtype] += 1

    total = sum(type_counts.values())

    with open(out_path, "w") as out:
        out.write("VARIANT_TYPE\tCOUNT\tFRACTION\n")
        for vtype in sorted(type_counts):
            count = type_counts[vtype]
            frac = count / total if total > 0 else 0
            out.write(f"{vtype}\t{count}\t{frac:.6f}\n")

        out.write(f"TOTAL_VARIANTS\t{total}\t1.0\n")

def main():
    ap = argparse.ArgumentParser(
        description="Summarize lead loci with unique/boosted LD and 'novel by override-only nominal' definition."
    )
    ap.add_argument("--ld", required=True, help="PLINK --r2 gz output (leadLD.ld.gz)")
    ap.add_argument("--leads-table", required=True, help="OUT.leads.from_clumped.tsv")
    ap.add_argument("--clump-table", required=True, help="OUT.for_clump.tsv (must contain MarkerID and P_CLUMP)")
    ap.add_argument("--clumped", required=True, help="OUT.clumped (to get members via SP2)")
    ap.add_argument("--boosted", required=True, help="Boosted IDs (one per line)")
    ap.add_argument("--unique", required=True, help="Unique IDs (one per line or table with MarkerID column)")
    ap.add_argument("--r2-link", type=float, default=0.8, help="Main LD threshold for detailed outputs (default 0.8)")
    ap.add_argument("--p-threshold", type=float, default=5e-5, help="Nominal p threshold (default 5e-5)")
    ap.add_argument("--out-prefix", required=True, help="Output prefix")
    ap.add_argument("--ld-sweep", action="store_true", help="Also compute LD sweep table")
    ap.add_argument("--ld-sweep-min", type=float, default=0.2, help="Sweep min (default 0.2)")
    ap.add_argument("--ld-sweep-max", type=float, default=1.0, help="Sweep max (default 1.0)")
    ap.add_argument("--ld-step", type=float, default=0.05, help="Sweep step (default 0.05)")
    args = ap.parse_args()

    out_prefix = args.out_prefix
    r2_link = float(args.r2_link)
    p_thr = float(args.p_threshold)

    # Load lists
    boosted = load_boosted_ids(args.boosted)
    unique = load_id_list(args.unique)
    override = unique | boosted

    # Leads
    lead_order, lead_info = parse_leads_table(args.leads_table)
    lead_set = set(lead_order)

    # Clumped members
    members_by_lead, all_members = parse_clumped_members(args.clumped)

    # P_CLUMP for members only
    p_clump = load_p_clump_for_members(args.clump_table, all_members)

    # nominal member counts
    nominal_counts = compute_nominal_member_counts(members_by_lead, p_clump, override, p_thr)

    # LD counts + max r2
    unique_hits, boosted_hits, maxU, maxB = parse_ld_counts_and_max(
        args.ld, lead_set, unique, boosted, r2_link
    )

    # Outputs
    locus_stats = f"{out_prefix}.locus_stats.r2ge{r2_link}.tsv"
    overall_stats = f"{out_prefix}.overall_stats.r2ge{r2_link}.tsv"

    write_locus_stats(
        locus_stats, lead_order, lead_info, unique, boosted,
        unique_hits, boosted_hits, nominal_counts
    )
    write_overall_stats(overall_stats, locus_stats, r2_link, p_thr)

    variant_type_out = f"{out_prefix}.variant_type_breakdown.r2ge{r2_link}.tsv"

    compute_variant_type_breakdown_from_ids(
        variant_type_out,
        lead_order,
        unique_hits,
        boosted_hits,
        unique,
        boosted,
    )

    print(f"  Variant type breakdown: {variant_type_out}")

    if args.ld_sweep:
        sweep_out = f"{out_prefix}.ld_sweep_counts.tsv"
        write_ld_sweep(
            sweep_out,
            lead_order,
            unique,
            boosted,
            maxU,
            maxB,
            float(args.ld_sweep_min),
            float(args.ld_sweep_max),
            float(args.ld_step),
        )

    print("Done.")
    print(f"  Locus stats:   {locus_stats}")
    print(f"  Overall stats: {overall_stats}")
    if args.ld_sweep:
        print(f"  LD sweep:      {out_prefix}.ld_sweep_counts.tsv")


if __name__ == "__main__":
    main()