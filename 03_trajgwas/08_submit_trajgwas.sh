#!/bin/bash
# Submission helper for the TrajGWAS stage (collects the commands we ran on the cluster).
# Run from the working directory $WHI_DIR/10_Rebuild_v2 after copying the scripts from 03_trajgwas/,
# 05_results_tables_figures/ and config.sh into it, and after 06_trajgwas_setup.sh has finished.
# Usage: bash 08_submit_trajgwas.sh full | lookup | leads | figures
#   full    : genome-wide score tests, 22 chromosomes x 5 group/site sets, then the summary job
#   lookup  : rs4233949 (chr2:54,659,707) score and Wald tests. This variant fails R2 >= 0.8 in PAGE, so it is
#             tested from the R2 >= 0.3 hard calls (WHI_merged_v3) restricted to the analysed women.
#   leads   : group-specific QC (MAF >= 0.01, MAC >= 20, HWE P >= 1e-6), lead loci, Wald tests for the leads
#   figures : Manhattan and QQ plots
set -euo pipefail
source "$PWD/config.sh"
SBATCH_R="sbatch --parsable --time=01:00:00 --mem=64G --account=YOUR_ACCOUNT"
R_ENV="module load gcc/12.3.0 R/4.4.0"   # edit for your cluster
SETS="black hipbmd $PHENO_DIR/01_Black/pheno_black_hipbmd_data.csv
black spinebmd $PHENO_DIR/01_Black/pheno_black_spinebmd_data.csv
hispanic hipbmd $PHENO_DIR/03_Hispanic/pheno_hispanic_hipbmd_data.csv
hispanic spinebmd $PHENO_DIR/03_Hispanic/pheno_hispanic_spinebmd_data.csv
ai_an hipbmd $PHENO_DIR/04_AI_AN/pheno_ai_hipbmd_data.csv"
# (The AI/AN spine set was attempted but no null model converged; see the Methods.)

case "${1:-}" in
full)
  : > full_jobs.txt
  while read g o pheno; do
    J=$(sbatch --parsable --array=1-22 07_trajgwas_run.sh $g $o $pheno); echo "$g $o $J" | tee -a full_jobs.txt
  done <<< "$SETS"
  sbatch --dependency=afterany:$(awk '{printf "%s%s", (NR>1?":":""), $3}' full_jobs.txt) 09_summary.sh ;;
lookup)
  mkdir -p lookup
  awk '{print $1, $2}' WHI_v5_r08.fam > lookup/keep_v5.txt
  plink --bfile WHI_merged_v3 --keep lookup/keep_v5.txt --chr 2 --from-bp 54659707 --to-bp 54659707 --make-bed --out lookup/sptbn1_rs4233949 > /dev/null
  awk '{$1=$2; print}' lookup/sptbn1_rs4233949.fam > lookup/tmp.fam && mv lookup/tmp.fam lookup/sptbn1_rs4233949.fam
  cut -f2 lookup/sptbn1_rs4233949.bim > lookup/sptbn1_snp.txt
  while read g o pheno; do
    GENOPREFIX=$PWD/lookup/sptbn1_rs4233949 sbatch --array=2 --job-name=traj_lookup 07_trajgwas_run.sh $g $o $pheno
    GENOPREFIX=$PWD/lookup/sptbn1_rs4233949 WALD=$PWD/lookup/sptbn1_snp.txt sbatch --array=2 --job-name=traj_lookup 07_trajgwas_run.sh $g $o $pheno
  done <<< "$SETS" ;;
leads)
  sbatch --wait --time=01:00:00 --mem=64G --account=YOUR_ACCOUNT --wrap="$R_ENV; Rscript 13_group_qc.R"
  while read g o pheno; do
    L=$PWD/trajgwas/summary/wald_leads_${g}_${o}.txt; [ -s $L ] || { echo "no lead list for $g $o"; continue; }
    A=$(cut -d: -f1 $L | sort -nu | paste -sd, -)
    WALD=$L sbatch --array=$A --job-name=traj_wald 07_trajgwas_run.sh $g $o $pheno
  done <<< "$SETS" ;;
figures)
  $SBATCH_R --wrap="$R_ENV; Rscript 12_figures.R" ;;
*) echo "usage: bash 08_submit_trajgwas.sh full | lookup | leads | figures"; exit 1 ;;
esac
