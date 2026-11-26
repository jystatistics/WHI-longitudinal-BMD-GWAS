#!/bin/bash
#SBATCH --job-name=dedup_v3
#SBATCH --output=dedup_v3_%j.out
#SBATCH --time=06:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --mail-type=END,FAIL
# 2026-11-25. WHI_merged_v2 has 3,561 rows but only 2,649 unique women: 912 women were genotyped in both
# PAGE and SHARE and appear twice (same IID, different FID). This job
#  1. confirms each same-IID pair is the same person (IBD PI_HAT on the pruned SNPs),
#  2. keeps the SHARE copy (genome-wide Affymetrix 6.0 array) and drops the PAGE copy -> WHI_merged_v3,
#  3. checks relatedness among the unique women, and redoes the PCA without duplicates -> WHI_PCA_v3.
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # site paths (see config.sh)
cd ${WHI_DIR}/10_Rebuild_v2
awk '{print $2}' WHI_merged_v2.fam | sort | uniq -d > dup_iid.txt
cat PAGE_c1_v2.fam PAGE_c2_v2.fam | awk 'NR==FNR{d[$1];next} ($2 in d){print $1,$2}' dup_iid.txt - > remove_page_dups.txt
echo "duplicated IIDs: $(wc -l < dup_iid.txt); PAGE rows to drop: $(wc -l < remove_page_dups.txt)"
[ $(wc -l < dup_iid.txt) -eq $(wc -l < remove_page_dups.txt) ] || { echo "Each duplicated IID should have exactly one PAGE row"; exit 1; }

# 1. same person? (PI_HAT ~1 for duplicates)
awk 'NR==FNR{d[$1];next} ($2 in d){print $1,$2}' dup_iid.txt WHI_merged_v2.fam > dup_rows.txt
plink --memory 120000 --threads 8 --bfile WHI_merged_v2 --keep dup_rows.txt --extract prune.prune.in --genome --min 0.5 --out dup_ibd > /dev/null
awk 'NR>1 && $2==$4 && $10>=0.9 {n++} END{print "same-IID pairs with PI_HAT>=0.9: " n+0}' dup_ibd.genome
awk 'NR>1 && $2==$4 {print $10}' dup_ibd.genome | sort -g | awk '{v[NR]=$1} END{if(NR) print "PI_HAT of same-IID pairs: min " v[1] ", median " v[int((NR+1)/2)] ", max " v[NR]}'

# 2. unique women
plink --memory 120000 --threads 8 --bfile WHI_merged_v2 --remove remove_page_dups.txt --make-bed --out WHI_merged_v3 > /dev/null
echo "WHI_merged_v3: $(wc -l < WHI_merged_v3.bim) variants, $(wc -l < WHI_merged_v3.fam) women, unique IIDs $(awk '{print $2}' WHI_merged_v3.fam | sort -u | wc -l)"
plink --memory 120000 --bfile WHI_merged_v3 --snp 2:54659707:C:G --freq --out sptbn1_v3 > /dev/null; cat sptbn1_v3.frq

# 3. relatedness among unique women, then PCA
plink --memory 120000 --threads 8 --bfile WHI_merged_v3 --maf 0.05 --geno 0.02 --hwe 1e-10 --indep-pairwise 200 50 0.2 --out prune_v3 > /dev/null
echo "pruned SNPs for PCA/IBD: $(wc -l < prune_v3.prune.in)"
plink --memory 120000 --threads 8 --bfile WHI_merged_v3 --extract prune_v3.prune.in --genome --min 0.1 --out rel_v3 > /dev/null
awk 'NR>1 {if($10>=0.9) d++; else if($10>=0.354) f++; else if($10>=0.177) s++} END{print "pairs: PI_HAT>=0.9 (duplicate) " d+0 ", 0.354-0.9 (1st degree) " f+0 ", 0.177-0.354 (2nd degree) " s+0}' rel_v3.genome
plink --memory 120000 --threads 8 --bfile WHI_merged_v3 --extract prune_v3.prune.in --pca 20 --out WHI_PCA_v3 > /dev/null
echo "WHI_PCA_v3.eigenvec: $(wc -l < WHI_PCA_v3.eigenvec) women"
echo "Done."
