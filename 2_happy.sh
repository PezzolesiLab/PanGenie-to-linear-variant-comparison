#!/usr/bin/env bash
source "${PIPELINE_ROOT}/lib/task_common.sh"

load_module_if_available bcftools
load_module_if_available miniconda3/latest

if command -v conda >/dev/null 2>&1; then
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate hap.py
fi

linear_truth="$(require_single_match 'final_combined_results_UPIDS_linear_*.vcf')"
pangenie_query="$(require_single_match 'final_combined_results_UPIDSpangenie_*.vcf')"
require_nonempty_file "$HAPPY_REF"
export HGREF="$HAPPY_HGREF"

log_msg "Sorting VCFs for hap.py"
bcftools sort "$pangenie_query" -Ov -o pangenie_sorted.vcf
bcftools sort "$linear_truth" -Ov -o linear_sorted.vcf

threads="${SLURM_CPUS_ON_NODE:-${SLURM_NTASKS:-1}}"
log_msg "Running hap.py with $threads threads"
python2.7 "$(which hap.py)" \
  linear_sorted.vcf \
  pangenie_sorted.vcf \
  --threads "$threads" \
  -r "$HAPPY_REF" \
  --window-size "$HAPPY_WINDOW_SIZE" \
  --write-vcf \
  --engine xcmp \
  --preprocess-truth \
  --no-roc \
  -o happy.Pangenie_VS_linear

log_msg "Done"
