#!/bin/bash
#SBATCH --job-name=dedup_v5
#SBATCH --output=dedup_v5_%j.out
#SBATCH --time=06:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --mail-type=END,FAIL
# 2026-11-25. Third duplicate pass. In WHI_v4_r08, 152 one-to-one pairs (114 SHARE-WHISP, 19 PAGE-WHISP,
# 12 GARNET-WHISP, 5 GARNET-PAGE, 1 SHARE-SHARE, 1 PAGE-PAGE) have IBD Z2 >= 0.5 with low Z1 (same woman,
# one noisy copy; relatives have Z2 <= ~0.25), and the Black phenotype file holds exactly one of the two IDs
# for 149 of them (never both). Drop the copy from the lower-priority batch (GARNET > WHIMS > SHARE > PAGE > WHISP)
# -> WHI_v5_r08; append to the ID alias list; redo PCA and the relatedness summary.
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # site paths (see config.sh)
cd ${WHI_DIR}/10_Rebuild_v2
awk 'NR==FNR{c[$1" "$2]=$3; next}
     FNR>1 && ($10>=0.9 || $9>=0.5) {
       p["GARNET"]=1; p["WHIMS"]=2; p["SHARE"]=3; p["PAGE"]=4; p["WHISP"]=5
       a=c[$1" "$2]; b=c[$3" "$4]
       if (p[a] <= p[b]) print $3, $4, $1, $2, b, a; else print $1, $2, $3, $4, a, b }' cohort_by_fidiid.txt rel_v4.genome > dup3_pairs.txt
awk '{print $1, $2}' dup3_pairs.txt | sort -u > remove_dup3.txt
echo "pairs: $(wc -l < dup3_pairs.txt); rows to drop: $(wc -l < remove_dup3.txt)"
awk '{print "drop " $5 " keep " $6}' dup3_pairs.txt | sort | uniq -c
cat alias_dup2.txt > alias_all.txt; awk '{print $2, $4}' dup3_pairs.txt >> alias_all.txt
echo "alias_all.txt: $(wc -l < alias_all.txt) dropped IDs"
plink --memory 120000 --threads 8 --bfile WHI_v4_r08 --remove remove_dup3.txt --make-bed --out WHI_v5_r08 > /dev/null
echo "WHI_v5_r08: $(wc -l < WHI_v5_r08.bim) variants, $(wc -l < WHI_v5_r08.fam) women"
awk 'NR==FNR{c[$1" "$2]=$3; next} {print c[$1" "$2]}' cohort_by_fidiid.txt WHI_v5_r08.fam | sort | uniq -c
plink --memory 120000 --threads 8 --bfile WHI_v5_r08 --maf 0.05 --geno 0.02 --hwe 1e-10 --exclude range longLD_hg19.txt --indep-pairwise 200 50 0.1 --out prune_v5 > /dev/null
echo "PCA SNPs: $(wc -l < prune_v5.prune.in)"
plink --memory 120000 --threads 8 --bfile WHI_v5_r08 --extract prune_v5.prune.in --pca 20 --out WHI_PCA_v5_r08 > /dev/null
awk '{for(j=3;j<=12;j++){s[j]+=$j; q[j]+=$j*$j}; n++; v[n]=$0} END{for(j=3;j<=12;j++){m[j]=s[j]/n; sd[j]=sqrt(q[j]/n-m[j]^2)}
  for(i=1;i<=n;i++){split(v[i],x," "); for(j=3;j<=12;j++) if((x[j]-m[j])/sd[j]>6||(x[j]-m[j])/sd[j]<-6) out[j-2]++}
  printf "PC outliers (|z|>6) by PC1-10:"; for(p=1;p<=10;p++) printf " %d", out[p]; print ""}' WHI_PCA_v5_r08.eigenvec
plink --memory 120000 --threads 8 --bfile WHI_v5_r08 --extract prune_v5.prune.in --genome --min 0.3 --out rel_v5 > /dev/null
awk 'NR==FNR{c[$1" "$2]=$3; next} FNR>1 {a=c[$1" "$2]; b=c[$3" "$4]; k=(a<b)?a"-"b:b"-"a; cl=($10>=0.9)?"dup":($10>=0.354)?"1st":"0.3-0.354"; x=cl" "k; n[x]++; z1[x]+=$8; z2[x]+=$9; if($9>=0.5) h[x]++}
  END{for(x in n) printf "%-28s pairs %5d  mean Z1 %.2f  mean Z2 %.2f  Z2>=0.5 %d\n", x, n[x], z1[x]/n[x], z2[x]/n[x], h[x]+0}' cohort_by_fidiid.txt rel_v5.genome | sort
echo "Done."
