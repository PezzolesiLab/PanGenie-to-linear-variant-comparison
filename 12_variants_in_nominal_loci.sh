#!/usr/bin/env bash
set -euo pipefail

# combined_clump_and_leadLD.sh
#
# PURPOSE (updated):
#   - Run PLINK clumping on a *combined* PLINK callset to get 1 lead variant per nominal locus.
#   - Compute lead↔nearby-variant LD (r^2) on the same combined callset.
#
# P-VALUE RULE:
#   - Use PanGenie GWAS p.value for every variant where PanGenie has the variant.
#   - Use Linear GWAS p.value only when the variant is absent/missing from PanGenie.
#
# This makes PanGenie the default p-value source for lead selection. Leads are simply whatever
# clumping selects given the P_CLUMP values defined above.
#
# INPUTS:
#   --bfile          Combined PLINK prefix (bed/bim/fam)
#   --pangenie-gwas  PanGenie GWAS results (MarkerID + p.value) [gz ok]
#   --linear-gwas    Linear GWAS results (MarkerID + p.value)   [gz ok]
#   --boosted-ids    one MarkerID per line (first field used)
#   --unique         PanGenie-unique file: either one MarkerID per line OR a table with "MarkerID" column
#
# OUTPUTS (prefix = OUT):
#   OUT.override.ids.txt               (union of unique + boosted IDs)
#   OUT.for_clump.tsv                  (CHR BP MarkerID P_CLUMP P_PANGENIE P_LINEAR SOURCE A1 A2)
#   OUT.clumped
#   OUT.leads.txt
#   OUT.leads.from_clumped.tsv
#   OUT.leads.from_bim.tsv
#   OUT.leads.pvals.tsv                (lead pvals: P_CLUMP/P_PANGENIE/P_LINEAR/SOURCE)
#   OUT.leadLD.ld.gz
#
# IMPORTANT:
#   All PLINK commands include --keep-allele-order.

usage() {
  cat <<'EOF'
Usage:
  combined_clump_and_leadLD.sh \
    --bfile COMBINED_PREFIX \
    --pangenie-gwas pangenie_gwas.tsv[.gz] \
    --linear-gwas   linear_gwas.tsv[.gz] \
    --boosted-ids   boosted.ids \
    --unique        pangenie_unique.{ids|tsv} \
    [--out OUT_PREFIX] \
    [--p-lead 0.05] \
    [--r2-indep 0.1] \
    [--kb 1000] \
    [--ld-min-r2 0] \
    [--ld-window 10000000] \
    [--threads N] \
    [--plink plink]

Required:
  --bfile          Combined PLINK prefix (bed/bim/fam)
  --pangenie-gwas  PanGenie GWAS file (MarkerID + p.value; gz ok)
  --linear-gwas    Linear GWAS file (MarkerID + p.value; gz ok)
  --boosted-ids    Boosted shared IDs (one per line; first field used)
  --unique         PanGenie-unique IDs file:
                    - either one MarkerID per line, OR
                    - a whitespace/tab table with header containing "MarkerID"

Optional:
  --out            Output prefix (default derived from inputs + params)
  --p-lead         Lead/index p-value threshold for clumping (default 0.05)
  --r2-indep       Lead-lead independence r^2 for clumping (default 0.1)
  --kb             Window size in kb for clumping/LD around lead (default 1000)
  --ld-min-r2      Minimum r^2 to report in LD table (default 0; can be huge output)
  --ld-window      Max variant-count window cap for LD (default 10000000)
  --threads        PLINK threads (default: SLURM_NTASKS if set, else 1)
  --plink          PLINK executable (default "plink")

Example:
  ./combined_clump_and_leadLD.sh \
    --bfile COMBINED.linear_plus_pangenie_unique_boosted \
    --pangenie-gwas pangenie.gwas.tsv.gz \
    --linear-gwas linear.gwas.tsv.gz \
    --boosted-ids 10_manhattan_plot.boosted_shared.pangenie-stronger.ALL.ids \
    --unique 4_compare_gwas_unique_buff10.only_in_pangenie.tsv \
    --p-lead 5e-5 --r2-indep 0.1 --kb 1000 --ld-min-r2 0.01 \
    --out phenotypeX.combined
EOF
}

die(){ echo "ERROR: $*" >&2; exit 2; }

BFILE=""
GWAS_PAN=""
GWAS_LIN=""
BOOSTED_FILE=""
UNIQUE_FILE=""
OUT=""
P_LEAD="0.05"
R2_INDEP="0.1"
KB="1000"
LD_MIN_R2="0"
LD_WINDOW="10000000"
THREADS="${SLURM_NTASKS:-1}"
PLINK_BIN="plink"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bfile)         BFILE="$2"; shift 2;;
    --pangenie-gwas) GWAS_PAN="$2"; shift 2;;
    --linear-gwas)   GWAS_LIN="$2"; shift 2;;
    --boosted-ids)   BOOSTED_FILE="$2"; shift 2;;
    --unique)        UNIQUE_FILE="$2"; shift 2;;
    --out)           OUT="$2"; shift 2;;
    --p-lead)        P_LEAD="$2"; shift 2;;
    --r2-indep)      R2_INDEP="$2"; shift 2;;
    --kb)            KB="$2"; shift 2;;
    --ld-min-r2)     LD_MIN_R2="$2"; shift 2;;
    --ld-window)     LD_WINDOW="$2"; shift 2;;
    --threads)       THREADS="$2"; shift 2;;
    --plink)         PLINK_BIN="$2"; shift 2;;
    -h|--help)       usage; exit 0;;
    *) die "Unknown arg: $1 (use --help)";;
  esac
done

[[ -n "$BFILE" ]]       || { usage; die "--bfile is required"; }
[[ -n "$GWAS_PAN" ]]    || { usage; die "--pangenie-gwas is required"; }
[[ -n "$GWAS_LIN" ]]    || { usage; die "--linear-gwas is required"; }
[[ -n "$BOOSTED_FILE" ]]|| { usage; die "--boosted-ids is required"; }
[[ -n "$UNIQUE_FILE" ]] || { usage; die "--unique is required"; }

[[ -f "${BFILE}.bed" ]] || die "Missing ${BFILE}.bed"
[[ -f "${BFILE}.bim" ]] || die "Missing ${BFILE}.bim"
[[ -f "${BFILE}.fam" ]] || die "Missing ${BFILE}.fam"
[[ -f "$GWAS_PAN" ]]    || die "Missing PanGenie GWAS file: $GWAS_PAN"
[[ -f "$GWAS_LIN" ]]    || die "Missing Linear GWAS file: $GWAS_LIN"
[[ -f "$BOOSTED_FILE" ]]|| die "Missing boosted IDs file: $BOOSTED_FILE"
[[ -f "$UNIQUE_FILE" ]] || die "Missing unique file: $UNIQUE_FILE"

command -v "$PLINK_BIN" >/dev/null 2>&1 || die "Cannot find PLINK executable: $PLINK_BIN"
command -v awk >/dev/null 2>&1 || die "Cannot find awk"
command -v sort >/dev/null 2>&1 || die "Cannot find sort"
command -v gzip >/dev/null 2>&1 || die "Cannot find gzip"
command -v mktemp >/dev/null 2>&1 || die "Cannot find mktemp"
command -v gunzip >/dev/null 2>&1 || die "Cannot find gunzip"

# Derive OUT if not provided
if [[ -z "$OUT" ]]; then
  basep="$(basename "$GWAS_PAN")"; basep="${basep%.gz}"; basep="${basep%.*}"
  basel="$(basename "$GWAS_LIN")"; basel="${basel%.gz}"; basel="${basel%.*}"
  OUT="combined.${basep}__${basel}.p${P_LEAD}.r2indep${R2_INDEP}.kb${KB}"
fi

tmp_pan=""
tmp_lin=""
cleanup() {
  local status=$?
  [[ -n "${tmp_pan}" && -f "${tmp_pan}" ]] && rm -f "${tmp_pan}"
  [[ -n "${tmp_lin}" && -f "${tmp_lin}" ]] && rm -f "${tmp_lin}"
  return "$status"
}
trap cleanup EXIT

PAN_IN="$GWAS_PAN"
LIN_IN="$GWAS_LIN"
if [[ "$GWAS_PAN" == *.gz ]]; then
  tmp_pan="$(mktemp "${OUT}.pangenie.XXXXXX.tsv")"
  gunzip -c "$GWAS_PAN" > "$tmp_pan"
  PAN_IN="$tmp_pan"
fi
if [[ "$GWAS_LIN" == *.gz ]]; then
  tmp_lin="$(mktemp "${OUT}.linear.XXXXXX.tsv")"
  gunzip -c "$GWAS_LIN" > "$tmp_lin"
  LIN_IN="$tmp_lin"
fi

# Output helper files
BOOSTED_IDS="${OUT}.boosted.ids.txt"
UNIQUE_IDS="${OUT}.unique.ids.txt"
OVERRIDE_IDS="${OUT}.override.ids.txt"
FOR_CLUMP="${OUT}.for_clump.tsv"

# -------------------------
# 0) Normalize boosted IDs (one per line; first field used) + dedupe
# -------------------------
awk 'NF>0 && $1!~/^#/ {gsub(/\r/,"",$1); gsub(/^"+|"+$/,"",$1); print $1}' "$BOOSTED_FILE" \
  | sort -u > "$BOOSTED_IDS"
[[ -s "$BOOSTED_IDS" ]] || die "Boosted IDs empty after parsing: $BOOSTED_IDS"

# -------------------------
# 1) Extract unique IDs
#    - if header has MarkerID => extract that column
#    - else treat as one ID per line
# -------------------------
first_nonempty="$(awk 'NF>0{print; exit}' "$UNIQUE_FILE")"

if echo "$first_nonempty" | awk 'BEGIN{FS="[ \t]+"} {for(i=1;i<=NF;i++) if($i=="MarkerID") found=1} END{exit(found?0:1)}'; then
  awk '
    BEGIN{FS="[ \t]+"; mid=0}
    NR==1{
      for(i=1;i<=NF;i++) if($i=="MarkerID") mid=i
      if(mid==0){ print "ERROR: MarkerID not found in header" > "/dev/stderr"; exit 2 }
      next
    }
    NF>0{
      id=$mid
      gsub(/\r/,"",id)
      gsub(/^"+|"+$/,"",id)
      if(id!="") print id
    }
  ' "$UNIQUE_FILE" | sort -u > "$UNIQUE_IDS"
else
  awk '
    NF==0{next}
    $0 ~ /^#/ {next}
    {id=$1; gsub(/\r/,"",id); gsub(/^"+|"+$/,"",id); if(id!="") print id}
  ' "$UNIQUE_FILE" | sort -u > "$UNIQUE_IDS"
fi

[[ -s "$UNIQUE_IDS" ]] || die "Unique IDs empty after parsing: $UNIQUE_IDS"

# -------------------------
# 2) Override set = unique ∪ boosted
# -------------------------
cat "$UNIQUE_IDS" "$BOOSTED_IDS" | sort -u > "$OVERRIDE_IDS"
[[ -s "$OVERRIDE_IDS" ]] || die "Override IDs empty: $OVERRIDE_IDS"

echo "[0/3] Building clumping input table from combined BIM + PanGenie-default p-values + Linear fallback..."

# -------------------------
# 3) Build FOR_CLUMP from:
#    - combined .bim (defines variant universe + chr/bp/alleles)
#    - PanGenie GWAS p.value whenever available
#    - Linear GWAS p.value only when PanGenie is missing
#
# Output:
#   CHR BP MarkerID P_CLUMP P_PANGENIE P_LINEAR SOURCE A1 A2
# SOURCE:
#   PAN_DEFAULT
#   PAN_DEFAULT_OVERRIDE_ID
#   PAN_MISSING_FALLBACK_LIN
#   PAN_MISSING_FALLBACK_LIN_OVERRIDE_ID
#   MISSING
# -------------------------
awk -v FS='[ \t]+' -v OFS="\t" -v ov="$OVERRIDE_IDS" '
  function norm(s){ gsub(/\r/,"",s); gsub(/^"+|"+$/,"",s); return s }
  function isnum(x){ return (x ~ /^[0-9.]+([eE][-+]?[0-9]+)?$/) }

  # --- File 1: combined BIM ---
  FILENAME==ARGV[1]{
    # .bim: CHR SNP CM BP A1 A2
    snp = norm($2)
    if(snp=="") next
    chr[snp]=$1
    bp[snp]=$4
    a1[snp]=$5
    a2[snp]=$6
    keep[snp]=1
    order[++n]=snp
    next
  }

  # --- Load override IDs (unique+boosted) ---
  BEGIN{
    while((getline line < ov) > 0){
      line = norm(line)
      if(line=="" || line ~ /^#/) continue
      override[line]=1
    }
    close(ov)
  }

  # --- File 2: PanGenie GWAS (store for any variant in combined BIM) ---
  FILENAME==ARGV[2]{
    if(FNR==1){
      mid=0; pv=0
      for(i=1;i<=NF;i++){ if($i=="MarkerID") mid=i; if($i=="p.value") pv=i }
      if(!mid || !pv){ print "ERROR: PanGenie GWAS missing MarkerID and/or p.value" > "/dev/stderr"; exit 2 }
      next
    }
    id = norm($(mid))
    if(!(id in keep)) next
    p = norm($(pv))
    if(p=="" || !isnum(p)) next
    if(p+0 <= 0) p=1e-300
    pPan[id]=p+0.0
    next
  }

  # --- File 3: Linear GWAS (store for any variant in combined BIM) ---
  FILENAME==ARGV[3]{
    if(FNR==1){
      mid=0; pv=0
      for(i=1;i<=NF;i++){ if($i=="MarkerID") mid=i; if($i=="p.value") pv=i }
      if(!mid || !pv){ print "ERROR: Linear GWAS missing MarkerID and/or p.value" > "/dev/stderr"; exit 2 }
      next
    }
    id = norm($(mid))
    if(!(id in keep)) next
    p = norm($(pv))
    if(p=="" || !isnum(p)) next
    if(p+0 <= 0) p=1e-300
    pLin[id]=p+0.0
    next
  }

  END{
    print "CHR","BP","MarkerID","P_CLUMP","P_PANGENIE","P_LINEAR","SOURCE","A1","A2"

    miss_pan=0; miss_lin=0; pan_default=0; fallback_lin=0; used_missing=0
    override_pan_default=0; override_fallback_lin=0

    for(i=1;i<=n;i++){
      id = order[i]
      if(!(id in keep)) continue

      pan = (id in pPan ? pPan[id] : "NA")
      lin = (id in pLin ? pLin[id] : "NA")
      if(pan == "NA") miss_pan++
      if(lin == "NA") miss_lin++

      p="NA"; src="MISSING"

      if(pan != "NA"){
        p = pan
        pan_default++
        if(id in override){
          src = "PAN_DEFAULT_OVERRIDE_ID"
          override_pan_default++
        } else {
          src = "PAN_DEFAULT"
        }
      } else if(lin != "NA"){
        p = lin
        fallback_lin++
        if(id in override){
          src = "PAN_MISSING_FALLBACK_LIN_OVERRIDE_ID"
          override_fallback_lin++
        } else {
          src = "PAN_MISSING_FALLBACK_LIN"
        }
      } else {
        used_missing++
      }

      # emit only if we have a p-value for clumping
      if(p=="NA") continue

      print chr[id], bp[id], id, p, pan, lin, src, a1[id], a2[id]
    }

    # Print diagnostics to stderr
    print "[INFO] FOR_CLUMP diagnostics:", \
          "override_ids=" length(override), \
          "bim_variants=" n, \
          "pan_missing=" miss_pan, \
          "lin_missing=" miss_lin, \
          "pangenie_default=" pan_default, \
          "fallback_to_linear=" fallback_lin, \
          "override_pangenie_default=" override_pan_default, \
          "override_fallback_to_linear=" override_fallback_lin, \
          "dropped_missing_both=" used_missing \
          > "/dev/stderr"
  }
' "${BFILE}.bim" "$PAN_IN" "$LIN_IN" > "$FOR_CLUMP"

[[ -s "$FOR_CLUMP" ]] || die "Clumping table empty: $FOR_CLUMP"

# Common PLINK opts (always keep allele order)
PLINK_OPTS=(--allow-extra-chr --keep-allele-order --threads "$THREADS")

echo "[1/3] Clumping to get lead variants per nominal locus..."
"$PLINK_BIN" "${PLINK_OPTS[@]}" \
  --bfile "$BFILE" \
  --clump "$FOR_CLUMP" \
  --clump-snp-field MarkerID \
  --clump-field "P_CLUMP" \
  --clump-p1 "$P_LEAD" \
  --clump-p2 1 \
  --clump-r2 "$R2_INDEP" \
  --clump-kb "$KB" \
  --out "$OUT"

[[ -s "${OUT}.clumped" ]] || die "Clumping produced no .clumped output: ${OUT}.clumped"

echo "[2/3] Extracting lead SNP list and lead tables..."

# Leads list (one SNP per locus)
awk '
  BEGIN{OFS="\t"}
  NR==1{next}
  $0 ~ /^#/ {next}
  NF>=5 {
    snp=$3
    if (!(snp in seen)) { seen[snp]=1; print snp }
  }
' "${OUT}.clumped" > "${OUT}.leads.txt"

[[ -s "${OUT}.leads.txt" ]] || die "No leads extracted (is --p-lead too strict?)."

# Lead locus table from .clumped:
# CHR  F  SNP  BP  P  TOTAL  NSIG ...
awk '
  BEGIN{OFS="\t"}
  NR==1{next}
  $0 ~ /^#/ {next}
  NF>=7 { print $1,$3,$4,$5,$6,$7 }
' "${OUT}.clumped" \
  | awk 'BEGIN{OFS="\t"} NR==1{print "CHR","SNP","BP","P_CLUMP","TOTAL","NSIG"} {print}' \
  > "${OUT}.leads.from_clumped.tsv"

# Lead BIM info from combined BIM
awk '
  BEGIN{OFS="\t"}
  NR==FNR {keep[$1]=1; next}
  ($2 in keep) {print $1,$4,$2,$5,$6}
' "${OUT}.leads.txt" "${BFILE}.bim" \
  | awk 'BEGIN{OFS="\t"} NR==1{print "CHR","BP","SNP","A1","A2"} {print}' \
  > "${OUT}.leads.from_bim.tsv"

# Lead p-values from FOR_CLUMP (P_CLUMP, P_PANGENIE, P_LINEAR, SOURCE)
awk -v FS='\t' -v OFS="\t" '
  NR==FNR{
    if(FNR==1) next
    pC[$3]=$4
    pP[$3]=$5
    pL[$3]=$6
    src[$3]=$7
    next
  }
  NR==1{ print "SNP","P_CLUMP","P_PANGENIE","P_LINEAR","SOURCE"; next }
  {
    snp=$1
    print snp, (snp in pC ? pC[snp] : "NA"),
               (snp in pP ? pP[snp] : "NA"),
               (snp in pL ? pL[snp] : "NA"),
               (snp in src ? src[snp] : "NA")
  }
' "$FOR_CLUMP" <(awk 'BEGIN{print "SNP"} {print $1}' "${OUT}.leads.txt") \
  > "${OUT}.leads.pvals.tsv"

echo "[3/3] Computing lead↔nearby-variant LD table (r^2) on the combined callset..."
"$PLINK_BIN" "${PLINK_OPTS[@]}" \
  --bfile "$BFILE" \
  --r2 gz \
  --ld-snp-list "${OUT}.leads.txt" \
  --ld-window "$LD_WINDOW" \
  --ld-window-kb "$KB" \
  --ld-window-r2 "$LD_MIN_R2" \
  --out "${OUT}.leadLD"

[[ -s "${OUT}.leadLD.ld.gz" ]] || die "LD output missing/empty: ${OUT}.leadLD.ld.gz"

echo
echo "Done."
echo "Override IDs (unique+boosted): ${OVERRIDE_IDS}"
echo "Clump input:                  ${FOR_CLUMP}"
echo "Leads:                        ${OUT}.leads.txt"
echo "Lead loci:                    ${OUT}.clumped"
echo "Lead table:                   ${OUT}.leads.from_clumped.tsv"
echo "Lead BIM info:                ${OUT}.leads.from_bim.tsv"
echo "Lead pvals:                   ${OUT}.leads.pvals.tsv"
echo "Lead LD pairs:                ${OUT}.leadLD.ld.gz"
