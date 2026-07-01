#!/usr/bin/env python3
import argparse
import gzip
import sys


def open_maybe_gzip(path, mode="rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def make_marker_id(chrom, pos, ref, alt):
    return f"{chrom}_{pos}_{ref}_{alt}"


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Extract shared, pangenie-unique, and linear-unique variant lists "
            "from a hap.py comparison VCF.\n\n"
            "Assumes:\n"
            "  TRUTH = linear presence VCF\n"
            "  QUERY = pangenie presence VCF\n"
        )
    )
    ap.add_argument("hap_vcf", help="hap.py output VCF (vcf or vcf.gz)")
    ap.add_argument("--prefix", default="hap", help="Output file prefix (default: hap)")

    args = ap.parse_args()

    shared = set()
    pangenie_unique = set()
    linear_unique = set()

    with open_maybe_gzip(args.hap_vcf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue

            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alt = fields[4]

            fmt = fields[8].split(":")
            truth = fields[9].split(":")
            query = fields[10].split(":")

            try:
                bd_idx = fmt.index("BD")
            except ValueError:
                sys.exit("ERROR: FORMAT field does not contain BD")

            bd_truth = truth[bd_idx]
            bd_query = query[bd_idx]

            mid = make_marker_id(chrom, pos, ref, alt)

            # Shared: TP on TRUTH side
            if bd_truth == "TP":
                shared.add(mid)

            # Linear-only: FN on TRUTH side
            elif bd_truth == "FN":
                linear_unique.add(mid)

            # Pangenie-only: FP on QUERY side
            elif bd_query == "FP":
                pangenie_unique.add(mid)

    # ---- Write outputs ----
    with open(f"{args.prefix}.shared_variants.txt", "w") as out:
        for v in sorted(shared):
            out.write(v + "\n")

    with open(f"{args.prefix}.pangenie_unique_variants.txt", "w") as out:
        for v in sorted(pangenie_unique):
            out.write(v + "\n")

    with open(f"{args.prefix}.linear_unique_variants.txt", "w") as out:
        for v in sorted(linear_unique):
            out.write(v + "\n")

    # ---- Summary ----
    print(f"Shared (TP):           {len(shared):,}")
    print(f"Pangenie-unique (FP):  {len(pangenie_unique):,}")
    print(f"Linear-unique (FN):    {len(linear_unique):,}")


if __name__ == "__main__":
    main()
