#!/usr/bin/env bash
set -euo pipefail

# build_combined_linear_plus_pangenie_unique_boosted.sh
#
# Creates a combined PLINK bed/bim/fam dataset:
#   COMBINED = (LINEAR minus boosted variants)  ∪  (PANGENIE subset to unique+boosted)
#
# IMPORTANT: All PLINK commands include --keep-allele-order.

usage() {
  cat <<'EOF'
Usage:
  build_combined_linear_plus_pangenie_unique_boosted.sh \
    --linear-prefix LINEAR \
    --pangenie-prefix PANGENIE \
    --boosted-ids 10_manhattan_plot.boosted_shared.pangenie-stronger.ALL.ids \
    --unique-table 4_compare_gwas_unique_buff10.only_in_pangenie.tsv \
    --out-prefix COMBINED.linear_plus_pangenie_unique_boosted \
    [--out-dir OUTDIR] \
    [--keep-intersect-samples] \
    [--plink plink] \
    [--debug]

Required:
  --linear-prefix     Prefix for linear callset PLINK files (no extension)
  --pangenie-prefix   Prefix for pangenie callset PLINK files (no extension)
  --boosted-ids       File with boosted IDs (one per line; first field used)
  --unique-table      Table containing a "MarkerID" column (unique variants)
  --out-prefix        Prefix name for combined output (no extension)

Optional:
  --out-dir OUTDIR             Output directory (default: .)
  --keep-intersect-samples     Intersect samples between callsets first
  --plink PLINK_BIN            PLINK executable name/path (default: plink)
  --debug                      Print extra diagnostics

EOF
}

die(){ echo "ERROR: $*" >&2; exit 2; }

LINEAR_PREFIX=""
PANGENIE_PREFIX=""
BOOSTED_FILE=""
UNIQUE_TABLE=""
OUT_PREFIX=""
OUT_DIR="."
KEEP_INTERSECT_SAMPLES=0
PLINK_BIN="plink"
DEBUG=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --linear-prefix) LINEAR_PREFIX="$2"; shift 2;;
    --pangenie-prefix) PANGENIE_PREFIX="$2"; shift 2;;
    --boosted-ids) BOOSTED_FILE="$2"; shift 2;;
    --unique-table) UNIQUE_TABLE="$2"; shift 2;;
    --out-prefix) OUT_PREFIX="$2"; shift 2;;
    --out-dir) OUT_DIR="$2"; shift 2;;
    --keep-intersect-samples) KEEP_INTERSECT_SAMPLES=1; shift 1;;
    --plink) PLINK_BIN="$2"; shift 2;;
    --debug) DEBUG=1; shift 1;;
    -h|--help) usage; exit 0;;
    *) die "Unknown arg: $1 (use --help)";;
  esac
done

[[ -n "$LINEAR_PREFIX" ]] || { usage; die "--linear-prefix is required"; }
[[ -n "$PANGENIE_PREFIX" ]] || { usage; die "--pangenie-prefix is required"; }
[[ -n "$BOOSTED_FILE" ]] || { usage; die "--boosted-ids is required"; }
[[ -n "$UNIQUE_TABLE" ]] || { usage; die "--unique-table is required"; }
[[ -n "$OUT_PREFIX" ]] || { usage; die "--out-prefix is required"; }

command -v "$PLINK_BIN" >/dev/null 2>&1 || die "Cannot find plink executable: $PLINK_BIN"
command -v awk >/dev/null 2>&1 || die "Cannot find awk"
command -v sort >/dev/null 2>&1 || die "Cannot find sort"
command -v comm >/dev/null 2>&1 || die "Cannot find comm"
command -v wc >/dev/null 2>&1 || die "Cannot find wc"
command -v mkdir >/dev/null 2>&1 || die "Cannot find mkdir"

mkdir -p "$OUT_DIR"

# Check PLINK files exist
for ext in bed bim fam; do
  [[ -f "${LINEAR_PREFIX}.${ext}" ]] || die "Missing: ${LINEAR_PREFIX}.${ext}"
  [[ -f "${PANGENIE_PREFIX}.${ext}" ]] || die "Missing: ${PANGENIE_PREFIX}.${ext}"
done
[[ -f "$BOOSTED_FILE" ]] || die "Missing boosted IDs file: $BOOSTED_FILE"
[[ -f "$UNIQUE_TABLE" ]] || die "Missing unique table: $UNIQUE_TABLE"

# Output paths (in OUT_DIR)
BOOSTED_IDS="${OUT_DIR}/boosted.ids"
UNIQUE_IDS="${OUT_DIR}/unique.ids"
PANGENIE_KEEP="${OUT_DIR}/pangenie.keep.unique_plus_boosted.ids"

LINEAR_WORK="${LINEAR_PREFIX}"
PANGENIE_WORK="${PANGENIE_PREFIX}"

PLINK_COMMON_OPTS=(--keep-allele-order)

# -------------------------
# Optional: intersect samples first
# -------------------------
if [[ "$KEEP_INTERSECT_SAMPLES" -eq 1 ]]; then
  if [[ $DEBUG -eq 1 ]]; then
    echo "[DEBUG] Intersecting samples between linear and pangenie..." >&2
  fi

  LINEAR_SAMPLES="${OUT_DIR}/linear.samples"
  PANGENIE_SAMPLES="${OUT_DIR}/pangenie.samples"
  SAMPLES_INTERSECT="${OUT_DIR}/samples.intersect"

  awk '{print $1"\t"$2}' "${LINEAR_PREFIX}.fam" | sort -k1,1 -k2,2 > "$LINEAR_SAMPLES"
  awk '{print $1"\t"$2}' "${PANGENIE_PREFIX}.fam" | sort -k1,1 -k2,2 > "$PANGENIE_SAMPLES"
  comm -12 "$LINEAR_SAMPLES" "$PANGENIE_SAMPLES" > "$SAMPLES_INTERSECT"

  nint=$(wc -l < "$SAMPLES_INTERSECT" | tr -d ' ')
  [[ "$nint" -gt 0 ]] || die "No intersecting samples found between callsets."

  if [[ $DEBUG -eq 1 ]]; then
    echo "[DEBUG] Intersect samples: ${nint}" >&2
  fi

  "$PLINK_BIN" "${PLINK_COMMON_OPTS[@]}" --bfile "$LINEAR_PREFIX" --keep "$SAMPLES_INTERSECT" --make-bed --out "${OUT_DIR}/LINEAR.common_samples" >/dev/null
  "$PLINK_BIN" "${PLINK_COMMON_OPTS[@]}" --bfile "$PANGENIE_PREFIX" --keep "$SAMPLES_INTERSECT" --make-bed --out "${OUT_DIR}/PANGENIE.common_samples" >/dev/null

  LINEAR_WORK="${OUT_DIR}/LINEAR.common_samples"
  PANGENIE_WORK="${OUT_DIR}/PANGENIE.common_samples"
fi

# -------------------------
# 0) Build boosted.ids (dedup)
# -------------------------
awk 'NF>0 && $1!~/^#/ {gsub(/\r/,"",$1); gsub(/^"+|"+$/,"",$1); print $1}' "$BOOSTED_FILE" \
  | sort -u > "$BOOSTED_IDS"
[[ -s "$BOOSTED_IDS" ]] || die "boosted.ids is empty after parsing: $BOOSTED_IDS"

# -------------------------
# 1) Extract unique.ids from unique table (MarkerID column)
# -------------------------
awk '
  BEGIN{FS="[ \t]+"; mid=0}
  NR==1{
    for(i=1;i<=NF;i++) if($i=="MarkerID") mid=i
    if(mid==0){print "ERROR: MarkerID column not found in unique table header" > "/dev/stderr"; exit 2}
    next
  }
  NF>0{
    id=$mid
    gsub(/\r/,"",id)
    gsub(/^"+|"+$/,"",id)
    if(id!="") print id
  }
' "$UNIQUE_TABLE" | sort -u > "$UNIQUE_IDS"
[[ -s "$UNIQUE_IDS" ]] || die "unique.ids is empty after parsing: $UNIQUE_IDS"

# -------------------------
# 2) Make PanGenie keep list = unique ∪ boosted
# -------------------------
cat "$UNIQUE_IDS" "$BOOSTED_IDS" | sort -u > "$PANGENIE_KEEP"
[[ -s "$PANGENIE_KEEP" ]] || die "PanGenie keep list is empty: $PANGENIE_KEEP"

if [[ $DEBUG -eq 1 ]]; then
  echo "[DEBUG] boosted.ids:  $(wc -l < "$BOOSTED_IDS")" >&2
  echo "[DEBUG] unique.ids:   $(wc -l < "$UNIQUE_IDS")" >&2
  echo "[DEBUG] keep.ids:     $(wc -l < "$PANGENIE_KEEP")" >&2
fi

# -------------------------
# 3) Subset PanGenie to keep list
# -------------------------
PANGENIE_SUB="${OUT_DIR}/PANGENIE.unique_plus_boosted"
"$PLINK_BIN" "${PLINK_COMMON_OPTS[@]}" --bfile "$PANGENIE_WORK" \
  --extract "$PANGENIE_KEEP" \
  --make-bed \
  --out "$PANGENIE_SUB" >/dev/null

# -------------------------
# 4) Remove boosted variants from linear
# -------------------------
LINEAR_SUB="${OUT_DIR}/LINEAR.no_boosted"
"$PLINK_BIN" "${PLINK_COMMON_OPTS[@]}" --bfile "$LINEAR_WORK" \
  --exclude "$BOOSTED_IDS" \
  --make-bed \
  --out "$LINEAR_SUB" >/dev/null

# -------------------------
# 5) Ensure no overlapping variant IDs remain (avoid merge conflicts)
# -------------------------
LINEAR_VIDS="${OUT_DIR}/linear.vids"
PANGENIE_VIDS="${OUT_DIR}/pangenie.vids"
OVERLAP_VIDS="${OUT_DIR}/overlap.vids"

cut -f2 "${LINEAR_SUB}.bim" | sort -u > "$LINEAR_VIDS"
cut -f2 "${PANGENIE_SUB}.bim" | sort -u > "$PANGENIE_VIDS"
comm -12 "$LINEAR_VIDS" "$PANGENIE_VIDS" > "$OVERLAP_VIDS"

nover=$(wc -l < "$OVERLAP_VIDS" | tr -d ' ')
if [[ $DEBUG -eq 1 ]]; then
  echo "[DEBUG] Overlapping variant IDs after exclusion: ${nover}" >&2
fi

LINEAR_MERGE_BASE="$LINEAR_SUB"
if [[ "$nover" -gt 0 ]]; then
  if [[ $DEBUG -eq 1 ]]; then
    echo "[DEBUG] Removing ${nover} overlapping IDs from linear to allow merge (PanGenie wins)..." >&2
  fi
  LINEAR_SUB2="${OUT_DIR}/LINEAR.no_boosted.no_overlap"
  "$PLINK_BIN" "${PLINK_COMMON_OPTS[@]}" --bfile "$LINEAR_SUB" \
    --exclude "$OVERLAP_VIDS" \
    --make-bed \
    --out "$LINEAR_SUB2" >/dev/null
  LINEAR_MERGE_BASE="$LINEAR_SUB2"
fi

# -------------------------
# 6) Merge into combined callset
# -------------------------
OUT_COMBINED="${OUT_DIR}/${OUT_PREFIX}"

"$PLINK_BIN" "${PLINK_COMMON_OPTS[@]}" --bfile "$LINEAR_MERGE_BASE" \
  --bmerge "${PANGENIE_SUB}.bed" "${PANGENIE_SUB}.bim" "${PANGENIE_SUB}.fam" \
  --make-bed \
  --out "$OUT_COMBINED" >/dev/null

# -------------------------
# Summary
# -------------------------
echo "Done."
echo "  boosted.ids:                     $BOOSTED_IDS ($(wc -l < "$BOOSTED_IDS"))"
echo "  unique.ids:                      $UNIQUE_IDS ($(wc -l < "$UNIQUE_IDS"))"
echo "  pangenie keep ids:               $PANGENIE_KEEP ($(wc -l < "$PANGENIE_KEEP"))"
echo "  PanGenie subset prefix:          $PANGENIE_SUB"
echo "  Linear subset prefix:            $LINEAR_SUB"
if [[ "$nover" -gt 0 ]]; then
  echo "  Overlap removed from linear:     $OVERLAP_VIDS (${nover})"
  echo "  Linear merge base prefix:        $LINEAR_MERGE_BASE"
fi
echo "  Combined output prefix:          $OUT_COMBINED"