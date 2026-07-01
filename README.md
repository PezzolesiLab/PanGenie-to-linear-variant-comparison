# PanGenie-to-linear-variant-comparison

This document summarizes the primary per-phenotype analysis workflow for the
PanGenie versus linear-reference comparison. It omits scheduler orchestration,
quality-control wrappers, plotting-only commands, aggregate summary commands,
logs, generated outputs, and local absolute filesystem paths.

Commands are shown as they are run from a phenotype-specific working directory.
Local inputs are represented with placeholders rather than machine-specific
paths.

## Placeholder Inputs

- `<analysis_scripts>`: directory containing the analysis scripts.
- `<pangenie_header_vcf>`: PanGenie VCF used to provide VCF header/reference metadata.
- `<linear_header_vcf>`: linear-reference VCF used to provide VCF header/reference metadata.
- `<reference_directory>`: directory containing reference files used by hap.py.
- `<concordance_regions_bed>`: BED file of high-concordance regions.
- `<annotation_beds>`: GIAB annotation BED files used for enrichment analyses.
- `<linear_plink_prefix>`: PLINK prefix for linear-reference genotypes.
- `<pangenie_plink_prefix>`: PLINK prefix for PanGenie genotypes.

## Common Tunable Parameter Values Used Below

| Parameter | Value |
| --- | --- |
| `HAPPY_WINDOW_SIZE` | `50` |
| `COMPARE_WINDOW_BP` | `10` |
| `MIN_CONCORDANCE` | `100` |
| `P_THRESHOLD` | `5e-5` |
| `P_FOLD` | `2` |
| `BOOST_DIRECTION` | `pangenie-stronger` |
| `P_MAX` | `1` |
| `QQ_MAX_POINTS` | `100000` |
| `CLUMP_P_LEAD` | `5e-5` |
| `CLUMP_R2_INDEP` | `0.1` |
| `CLUMP_KB` | `250` |
| `CLUMP_LD_MIN_R2` | `0.2` |
| `CLUMP_LD_WINDOW` | `10000000` |
| `LD_R2_LINK` | `0.8` |
| `LD_SWEEP_MIN` | `0.2` |
| `LD_SWEEP_MAX` | `1.0` |
| `LD_STEP` | `0.05` |

## 1. Prepare GWAS Results for hap.py

PanGenie and linear GWAS marker IDs are canonicalized, then each GWAS result
file is converted to a hap.py-compatible VCF.

```bash
python <analysis_scripts>/canonicalize_gwas_marker_ids.py \
  "$PANGENIE_GWAS" \
  "${PANGENIE_GWAS}.canonicalized" \
  --report "${PANGENIE_GWAS}.marker_id_case_collisions.tsv"

mv "${PANGENIE_GWAS}.canonicalized" "$PANGENIE_GWAS"

python <analysis_scripts>/canonicalize_gwas_marker_ids.py \
  "$LINEAR_GWAS" \
  "${LINEAR_GWAS}.canonicalized" \
  --report "${LINEAR_GWAS}.marker_id_case_collisions.tsv"

mv "${LINEAR_GWAS}.canonicalized" "$LINEAR_GWAS"

python <analysis_scripts>/1_Create_hap.py_vcf_from_GWAS_results.py \
  "$PANGENIE_GWAS" \
  --header-vcf <pangenie_header_vcf> \
  --sidecar-prefix 1_create_happy_pangenie \
  --remove-marker-ids "${PANGENIE_GWAS}.marker_id_case_collisions.tsv" \
  --sample-name FAKE \
  > "$PANGENIE_HAPPY_VCF"

python <analysis_scripts>/1_Create_hap.py_vcf_from_GWAS_results.py \
  "$LINEAR_GWAS" \
  --header-vcf <linear_header_vcf> \
  --sidecar-prefix 1_create_happy_linear \
  --remove-marker-ids "${LINEAR_GWAS}.marker_id_case_collisions.tsv" \
  --sample-name FAKE \
  > "$LINEAR_HAPPY_VCF"
```

## 2. Compare PanGenie and Linear VCFs with hap.py

The two hap.py input VCFs are sorted and compared using the linear-reference VCF
as truth and the PanGenie VCF as query. The pipeline wrapper containing this
external hap.py invocation is included in this bundle as `scripts/02_happy.sh`;
`scripts/1_Create_hap.py_vcf_from_GWAS_results.py` only prepares the input VCFs.

```bash
bcftools sort "$PANGENIE_HAPPY_VCF" -Ov -o pangenie_sorted.vcf
bcftools sort "$LINEAR_HAPPY_VCF" -Ov -o linear_sorted.vcf

export HGREF=<reference_directory>

python2.7 "$(which hap.py)" \
  linear_sorted.vcf \
  pangenie_sorted.vcf \
  --threads "$THREADS" \
  -r hs1.fa \
  --window-size 50 \
  --write-vcf \
  --engine xcmp \
  --preprocess-truth \
  --no-roc \
  -o happy.Pangenie_VS_linear
```

## 3. Extract Shared and Method-Unique Variant Sets

Variant sets are extracted from the hap.py comparison VCF.

```bash
python <analysis_scripts>/3_hap_py_extract_variant_sets.py \
  happy.Pangenie_VS_linear.vcf.gz \
  --prefix combined_run
```

## 4. Recover GWAS Rows for Shared and Unique Variants

PanGenie-unique, linear-unique, and shared variant IDs are linked back to the
original GWAS result rows.

```bash
python <analysis_scripts>/4_compare_gwas.py \
  "$PANGENIE_GWAS" \
  "$LINEAR_GWAS" \
  --pangenie-unique combined_run.pangenie_unique_variants.txt \
  --linear-unique combined_run.linear_unique_variants.txt \
  --shared-pre combined_run.shared_variants.txt \
  --prefix 4_compare_gwas_unique_buff10 \
  --window 10 \
  --pangenie-skipped-sv-ids 1_create_happy_pangenie.skipped_sv_ids.txt \
  --pangenie-removed-multisite-ids 1_create_happy_pangenie.removed_multisite_ids.txt \
  --linear-removed-multisite-ids 1_create_happy_linear.removed_multisite_ids.txt
```

## 5. Compare Shared-Variant Effect Sizes and P-values

Shared variants are filtered to high-concordance regions and compared for
PanGenie versus linear effect and association strength.

```bash
python <analysis_scripts>/5_beta_and_pvalue_correlation.py \
  --shared-pre-metrics 4_compare_gwas_unique_buff10.shared_variants_metrics.pre_concordance.tsv \
  --prefix 5_beta_and_pvalue_buff10 \
  --concordance-bed <concordance_regions_bed> \
  --min-concordance 100
```

## 6. Test Annotation Enrichment of PanGenie-Unique Variants

PanGenie-unique and linear-unique variants are compared across annotation
classes.

```bash
python <analysis_scripts>/6_unique_pangenie_enrichment.py \
  --pangenie-unique 4_compare_gwas_unique_buff10.only_in_pangenie.tsv \
  --linear-unique 4_compare_gwas_unique_buff10.only_in_linear.tsv \
  --shared-pre 4_compare_gwas_unique_buff10.shared_gwas_ids.pre_concordance.txt \
  --beds <annotation_beds> \
  --prefix 6_unique_pangenie
```

## 7. Test Annotation Enrichment of Shared P-value Outliers

Shared variants with stronger association evidence in one representation are
tested for enrichment across annotation classes.

```bash
python <analysis_scripts>/7_common_enrichment.py \
  "$PANGENIE_GWAS" \
  "$LINEAR_GWAS" \
  --pangenie-unique 4_compare_gwas_unique_buff10.only_in_pangenie.tsv \
  --linear-unique 4_compare_gwas_unique_buff10.only_in_linear.tsv \
  --shared-ids 5_beta_and_pvalue_buff10.shared_gwas_ids.txt \
  --beds <annotation_beds> \
  --p-fold 2 \
  --direction both \
  --prefix 7_common_pfold2
```

## 8. Define Boosted Shared Variants

Shared variants are filtered for PanGenie-stronger nominal associations and creates Manhattan/QQ plot files.

```bash
python <analysis_scripts>/10_plot_manhattan_linear_plus_pangenie.py \
  "$LINEAR_GWAS" \
  "$PANGENIE_GWAS" \
  --shared-ids 5_beta_and_pvalue_buff10.shared_gwas_ids.txt \
  --pangenie-unique-ids 4_compare_gwas_unique_buff10.only_in_pangenie.tsv \
  --prefix 10_manhattan_plot \
  --p-max 1 \
  --nominal-p 5e-5 \
  --p-fold 2 \
  --direction pangenie-stronger \
  --unique-point-size 20 \
  --boosted-point-size 20 \
  --nominal-tsv 10_manhattan_plot_nominal_variant_list.tsv \
  --qq-max-points 100000
```

## 9. Build a Combined PLINK Dataset for Locus-Level Analyses

The linear genotype set is combined with PanGenie variants that are unique or
nominally boosted. The boosted ID file represents shared variants meeting the
nominal p-value, fold-change, and direction criteria above.

```bash
bash <analysis_scripts>/12_create_union_plink_set.sh \
  --linear-prefix <linear_plink_prefix> \
  --pangenie-prefix <pangenie_plink_prefix> \
  --boosted-ids 10_manhattan_plot.boosted_shared.pangenie-stronger.ALL.ids \
  --unique-table 4_compare_gwas_unique_buff10.only_in_pangenie.tsv \
  --out-prefix 12_COMBINED.linear_plus_pangenie_unique_boosted \
  --out-dir 12_combined_linear_and_pangneieUniqueBoosted_plink \
  --keep-intersect-samples \
  --debug
```

## 10. Clump Nominal Loci and Compute LD

Nominal lead variants are clumped, and LD between lead variants and candidate
PanGenie unique/boosted variants is computed.

```bash
bash <analysis_scripts>/12_variants_in_nominal_loci.sh \
  --bfile 12_combined_linear_and_pangneieUniqueBoosted_plink/12_COMBINED.linear_plus_pangenie_unique_boosted \
  --pangenie-gwas "$PANGENIE_GWAS" \
  --linear-gwas "$LINEAR_GWAS" \
  --boosted-ids 10_manhattan_plot.boosted_shared.pangenie-stronger.ALL.ids \
  --unique 4_compare_gwas_unique_buff10.only_in_pangenie.tsv \
  --out 12_variants_in_nominal_loci \
  --p-lead 5e-5 \
  --r2-indep 0.1 \
  --kb 250 \
  --ld-min-r2 0.2 \
  --ld-window 10000000 \
  --threads "$THREADS"
```

## 11. Identify PanGenie Unique/Boosted Variants Linked to Lead Loci

Candidate PanGenie unique/boosted variants are evaluated for LD overlap with
lead loci, including an LD-threshold sweep.

```bash
python <analysis_scripts>/13_overlap.py \
  --ld 12_variants_in_nominal_loci.leadLD.ld.gz \
  --leads-table 12_variants_in_nominal_loci.leads.from_clumped.tsv \
  --clump-table 12_variants_in_nominal_loci.for_clump.tsv \
  --clumped 12_variants_in_nominal_loci.clumped \
  --boosted 10_manhattan_plot.boosted_shared.pangenie-stronger.ALL.ids \
  --unique 4_compare_gwas_unique_buff10.only_in_pangenie.tsv \
  --r2-link 0.8 \
  --p-threshold 5e-5 \
  --out-prefix 13_LD_overlap \
  --ld-sweep \
  --ld-sweep-min 0.2 \
  --ld-sweep-max 1.0 \
  --ld-step 0.05
```
