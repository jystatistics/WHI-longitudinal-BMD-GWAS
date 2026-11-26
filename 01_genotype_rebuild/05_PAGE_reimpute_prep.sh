#!/bin/bash
#SBATCH --job-name=page_reimp_prep
#SBATCH --output=page_reimp_prep_%j.out
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --mail-type=END,FAIL
# 2026-11-24. Michigan dropped whole 20 Mb chunks for PAGE because one sample had call rate <50% there
# (chunks-excluded.txt, "Low Sample Call Rates"). This finds those samples in each excluded chunk,
# removes them, and writes new upload VCFs for the affected chromosomes only.
# Output: 02_PAGE/05_Rebuild/reimpute/c1 and c2 (upload the PAGE_cX-chrN.vcf.gz files from there).
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # site paths (see config.sh)
P=${WHI_DBGAP_DIR}/02_PAGE
for CG in c1 c2; do
  S=$P/phg000910.v1.PAGE_WHI.genotype-calls-matrixfmt.$CG/filtered_subject_level_set
  X=$S/02_after_imputation/chunks-excluded.txt
  O=$P/05_Rebuild/reimpute/$CG; mkdir -p $O; cd $O
  echo "===== PAGE $CG"
  : > low_callrate.txt
  for ch in $(awk 'NR>1 && $1 ~ /^chunk_[0-9]+_/{print $1}' $X); do   # autosomes only (chrX chunks are not used)
    IFS=_ read -r _ c from to <<< "$ch"; from=$((10#$from)); to=$((10#$to))
    if plink --vcf $S/PAGE_${CG}-chr$c.vcf --double-id --chr $c --from-bp $from --to-bp $to --missing --out m_$ch > /dev/null 2>&1; then
      awk -v ch=$ch 'NR>1 && $6>0.5 {print $1, $2, ch, $6}' m_$ch.imiss >> low_callrate.txt
      echo "$ch: $(awk 'NR>1 && $6>0.5' m_$ch.imiss | wc -l) samples with call rate <50%; worst sample missing rate $(awk 'NR>1 && $6>m{m=$6} END{print m+0}' m_$ch.imiss)"
    else
      echo "$ch: plink could not read this region"
    fi
  done
  awk '{print $1, $2}' low_callrate.txt | sort -u > remove.txt
  echo "$CG: samples to remove: $(wc -l < remove.txt); of these in the 1,389 analysis women: $(awk 'NR==FNR{k[$2];next} ($2 in k)' $P/04_Merged/merged_page_cleaned.fam remove.txt | wc -l)"
  [ -s remove.txt ] || { echo "No low-call-rate samples found: stop and check before re-imputing"; continue; }
  for c in $(awk 'NR>1 && $1 ~ /^chunk_[0-9]+_/{split($1,a,"_"); print a[2]}' $X | sort -un); do
    plink --vcf $S/PAGE_${CG}-chr$c.vcf --double-id --keep-allele-order --remove remove.txt --recode vcf-iid bgz --out PAGE_${CG}-chr$c > /dev/null
    echo "chr$c: wrote PAGE_${CG}-chr$c.vcf.gz, $(zcat PAGE_${CG}-chr$c.vcf.gz | grep -vc '^#') variants, $(zcat PAGE_${CG}-chr$c.vcf.gz | grep -m1 '^#CHROM' | awk '{print NF-9}') samples"
  done
done
echo Done.
