#!/bin/bash
#SBATCH --job-name=fhs_prep
#SBATCH --output=fhs_prep_%A_%a.out
#SBATCH --time=06:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --mail-type=FAIL
#SBATCH --array=1-3
# 2026-11-25. FHS look-up of the rebuilt WHI findings (women only, same TrajGWAS model as WHI).
# Task 1 femoral neck, 2 trochanter, 3 lumbar spine L2-4 (FHS has no total hip).
# Steps: women with outcome -> LD-pruned SNPs -> drop relatives -> 10 PCs on unrelated women -> extract look-up SNPs.
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # site paths (see config.sh)
module load gcc/12.3.0 R/4.4.0
F=${FHS_DIR}; B=${WHI_DIR}/10_Rebuild_v2
O=$F/10_Lookup_v2; mkdir -p $O; cd $B
case $SLURM_ARRAY_TASK_ID in
  1) o=fnbmd; d=01_FNBMD; y=femoral_neck_bmd;;
  2) o=trbmd; d=02_TRBMD; y=trochanteric_bmd;;
  3) o=lsbmd; d=03_LSBMD; y=lumbar_spine_l24_bmd;;
esac
G=$F/$d/FHS_Merged_Female; P=$F/$d/fhs_pheno_female.csv
for f in $G.bed $G.bim $G.fam $P; do [ -f $f ] || { echo "MISSING $f"; exit 1; }; done
echo "== $o ($y): genotype $(wc -l < $G.fam) women, $(wc -l < $G.bim) variants; phenotype $P"
Rscript 16b_fhs_women.R $P $y $G.fam $O/women_$o.txt
PL="plink --threads 8 --memory 60000 --bfile $G"
$PL --keep $O/women_$o.txt --maf 0.05 --geno 0.02 --hwe 1e-6 --snps-only just-acgt --exclude range $B/longLD_hg19.txt \
    --indep-pairwise 200 50 0.1 --out $O/prune_$o > /dev/null
echo "pruned SNPs: $(wc -l < $O/prune_$o.prune.in)"
$PL --keep $O/women_$o.txt --extract $O/prune_$o.prune.in --genome --min 0.1 --out $O/rel_$o > /dev/null || true
Rscript 16d_fhs_unrelated.R $O/women_$o.txt $O/rel_$o.genome $O/unrel_$o.txt
$PL --keep $O/unrel_$o.txt --extract $O/prune_$o.prune.in --pca 10 --out $O/pca_$o > /dev/null
echo "PCs: $(wc -l < $O/pca_$o.eigenvec) women; eigenvalues: $(head -10 $O/pca_$o.eigenval | tr '\n' ' ')"
# look-up SNPs: rs4233949 (SPTBN1) + lead SNPs of the rebuilt WHI Black hip and spine scans (summary/leads_*.csv)
R2=$O/lookup_ranges.txt
{ echo "2 54659707 54659707 rs4233949"
  for s in black_hipbmd black_spinebmd; do L=$B/trajgwas/summary/leads_$s.csv
    [ -f $L ] && awk -F, -v s=$s 'NR>1 {print $3, $4, $4, s"_"NR-1}' $L || echo "WARNING no $L" >&2; done; } | sort -u -k1,1n -k2,2n > $R2.tmp.$o
mv $R2.tmp.$o $R2.$o
echo "look-up positions: $(wc -l < $R2.$o)"
$PL --keep $O/unrel_$o.txt --extract range $R2.$o --make-bed --out $O/lookup_$o > /dev/null
awk '{$1=$2; print}' $O/lookup_$o.fam > $O/lookup_$o.fam.tmp && mv $O/lookup_$o.fam.tmp $O/lookup_$o.fam   # FID = IID for TrajGWAS
plink --bfile $O/lookup_$o --freq --hardy --out $O/lookup_$o > /dev/null
echo "look-up genotype set: $(wc -l < $O/lookup_$o.fam) women, $(wc -l < $O/lookup_$o.bim) variants"
awk '$1==2 && $4==54659707' $O/lookup_$o.bim
