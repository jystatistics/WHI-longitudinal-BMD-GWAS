#!/bin/bash
#SBATCH --job-name=dedup_v4
#SBATCH --output=dedup_v4_%j.out
#SBATCH --time=06:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --mail-type=END,FAIL
# 2026-11-25. Second duplicate pass. On the R2>=0.8 SNPs, 134 pairs among the 2,649 women have PI_HAT>=0.9 and
# IBS>=0.977: 119 GARNET-WHISP and 15 GARNET-PAGE, i.e. the same woman under two different IDs.
# Keep the copy from the better-imputed batch (GARNET > WHIMS > SHARE > PAGE > WHISP; e.g. SPTBN1 R2 0.99 in
# GARNET vs 0.38 in WHISP), drop the other -> WHI_v4_r08. Writes alias_dup2.txt (dropped IID -> kept IID) so
# phenotype rows filed under the dropped ID still reach the kept genotypes; checks the phenotype files for women
# who appear under both IDs. Then PCA + relatedness on the deduplicated set.
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # site paths (see config.sh)
B=${WHI_DIR}
cd $B/10_Rebuild_v2
for f in GARNET_v2 WHIMS_v2 WHISP_v2 PAGE_c1_v2 PAGE_c2_v2 SHARE_c1_v2 SHARE_c2_v2; do awk -v c=${f%%_*} '{print $1, $2, c}' $f.fam; done > cohort_by_fidiid.txt
awk 'NR==FNR{c[$1" "$2]=$3; next}
     FNR>1 && $10>=0.9 && $12>=0.97 {
       p["GARNET"]=1; p["WHIMS"]=2; p["SHARE"]=3; p["PAGE"]=4; p["WHISP"]=5
       a=c[$1" "$2]; b=c[$3" "$4]
       if (p[a] <= p[b]) print $3, $4, $1, $2, b, a; else print $1, $2, $3, $4, a, b }' cohort_by_fidiid.txt rel_r08.genome > dup2_pairs.txt   # dropFID dropIID keepFID keepIID dropCohort keepCohort
awk '{print $1, $2}' dup2_pairs.txt | sort -u > remove_dup2.txt
awk '{print $2, $4}' dup2_pairs.txt > alias_dup2.txt
echo "duplicate pairs: $(wc -l < dup2_pairs.txt); rows to drop: $(wc -l < remove_dup2.txt)"
awk '{print "drop " $5 " keep " $6}' dup2_pairs.txt | sort | uniq -c
[ $(wc -l < remove_dup2.txt) -eq $(wc -l < dup2_pairs.txt) ] || echo "NOTE: some women are in more than one pair (check dup2_pairs.txt)"

# phenotype files: is a duplicated woman listed under one ID or both?
for g in 01_Black/pheno_black 03_Hispanic/pheno_hispanic 04_AI_AN/pheno_ai; do for o in hipbmd spinebmd; do
  f=$B/${g}_${o}_data.csv; [ -s $f ] || continue
  awk -F, -v o=$o 'NR==FNR{split($0,a," "); dm[a[1]]=a[2]; ks[a[2]]; next}
    FNR==1{for(i=1;i<=NF;i++){h=$i; gsub(/"/,"",h); if(h=="IID") ci=i; if(h=="time") ct=i; if(h==o) cy=i}; next}
    {id=$ci; gsub(/"/,"",id); r=$ct":"$cy
     if(id in dm){dk[dm[id]]=1; rd[dm[id]]=rd[dm[id]]" "r}
     if(id in ks){kk[id]=1; rk[id]=rk[id]" "r}}
    END{nb=0; ni=0; nd=0; nk=0; for(x in ks){ if((x in dk)&&(x in kk)){nb++; if(rd[x]==rk[x]) ni++} else if(x in dk) nd++; else if(x in kk) nk++ }
        printf "%s: duplicated women under kept ID only %d, dropped ID only %d, both IDs %d (identical records %d)\n", FILENAME, nk, nd, nb, ni}' \
    alias_dup2.txt $f
done; done

plink --memory 120000 --threads 8 --bfile WHI_v3_r08 --remove remove_dup2.txt --make-bed --out WHI_v4_r08 > /dev/null
echo "WHI_v4_r08: $(wc -l < WHI_v4_r08.bim) variants, $(wc -l < WHI_v4_r08.fam) women"
awk 'NR==FNR{c[$1" "$2]=$3; next} {print c[$1" "$2]}' cohort_by_fidiid.txt WHI_v4_r08.fam | sort | uniq -c

plink --memory 120000 --threads 8 --bfile WHI_v4_r08 --maf 0.05 --geno 0.02 --hwe 1e-10 --exclude range longLD_hg19.txt --indep-pairwise 200 50 0.1 --out prune_v4 > /dev/null
echo "PCA SNPs: $(wc -l < prune_v4.prune.in)"
plink --memory 120000 --threads 8 --bfile WHI_v4_r08 --extract prune_v4.prune.in --pca 20 --out WHI_PCA_v4_r08 > /dev/null
awk '{for(j=3;j<=12;j++){s[j]+=$j; q[j]+=$j*$j}; n++; v[n]=$0} END{for(j=3;j<=12;j++){m[j]=s[j]/n; sd[j]=sqrt(q[j]/n-m[j]^2)}
  for(i=1;i<=n;i++){split(v[i],x," "); for(j=3;j<=12;j++) if((x[j]-m[j])/sd[j]>6||(x[j]-m[j])/sd[j]<-6) out[j-2]++}
  printf "PC outliers (|z|>6) by PC1-10:"; for(p=1;p<=10;p++) printf " %d", out[p]; print ""}' WHI_PCA_v4_r08.eigenvec
plink --memory 120000 --threads 8 --bfile WHI_v4_r08 --extract prune_v4.prune.in --genome --min 0.3 --out rel_v4 > /dev/null
awk 'NR==FNR{c[$1" "$2]=$3; next} FNR>1 {a=c[$1" "$2]; b=c[$3" "$4]; k=(a<b)?a"-"b:b"-"a; cl=($10>=0.9)?"dup":($10>=0.354)?"1st":"0.3-0.354"; n[cl" "k]++; s[cl" "k]+=$12}
  END{for(x in n) printf "PI_HAT class/cohorts: %-28s pairs %5d  mean IBS %.3f\n", x, n[x], s[x]/n[x]}' cohort_by_fidiid.txt rel_v4.genome | sort
echo "Done."
