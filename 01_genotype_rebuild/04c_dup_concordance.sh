#!/bin/bash
#SBATCH --job-name=dup_conc
#SBATCH --output=dup_conc_%j.out
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --mail-type=END,FAIL
# 2026-11-25. The 912 women genotyped in both PAGE and SHARE agree on only ~92% of pruned genotypes (IBS).
# Compare each woman's PAGE copy with her SHARE copy call by call, and see WHERE they disagree:
# palindromic (A/T, C/G) SNPs, particular chromosomes (e.g. the re-imputed ones), MAF, or particular batches.
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # site paths (see config.sh)
cd ${WHI_DIR}/10_Rebuild_v2
cat SHARE_c1_v2.fam SHARE_c2_v2.fam | awk 'NR==FNR{d[$1];next} ($2 in d){print $1,$2}' dup_iid.txt - > share_dup_rows.txt
awk 'NR==FNR{f[$2]=$1;next} ($2 in f){print $1,$2,f[$2],$2}' share_dup_rows.txt remove_page_dups.txt > page_to_share_ids.txt
echo "pairs: SHARE copies $(wc -l < share_dup_rows.txt), PAGE copies mapped $(wc -l < page_to_share_ids.txt)"
# batch of each copy: IID, P/S, batch
{ for b in PAGE_c1 PAGE_c2; do awk -v b=$b 'NR==FNR{r[$1" "$2];next} (($1" "$2) in r){print $2, "P", b}' remove_page_dups.txt ${b}_v2.fam; done
  for b in SHARE_c1 SHARE_c2; do awk -v b=$b 'NR==FNR{r[$1" "$2];next} (($1" "$2) in r){print $2, "S", b}' share_dup_rows.txt ${b}_v2.fam; done; } > dup_copy_batch.txt
echo "copy batches: $(cut -d' ' -f3 dup_copy_batch.txt | sort | uniq -c | tr '\n' ' ')"

conc () {  # $1 = plink variant filter, $2 = label
  plink --memory 60000 --bfile WHI_merged_v2 --keep share_dup_rows.txt $1 --make-bed --out cc_share_$2 > /dev/null
  plink --memory 60000 --bfile WHI_merged_v2 --keep remove_page_dups.txt $1 --make-bed --out cc_page0_$2 > /dev/null
  plink --memory 60000 --bfile cc_page0_$2 --update-ids page_to_share_ids.txt --make-bed --out cc_page_$2 > /dev/null
  echo "===== $2: $(wc -l < cc_share_$2.bim) variants"
  plink --memory 60000 --bfile cc_share_$2 --bmerge cc_page_$2 --merge-mode 6 --out cc_$2 | grep -E 'overlapping|concordan' || true
  plink --memory 60000 --bfile cc_share_$2 --freq --out cc_share_$2 > /dev/null
  # per-SNP discordant calls (of ~912 pairs)
  awk 'NR>1{n[$1]++} END{for(s in n) print s, n[s]}' cc_$2.diff > cc_${2}_snp.txt
  awk -v N=$(wc -l < share_dup_rows.txt) 'FILENAME==ARGV[1]{d[$1]=$2; next} FNR==1{next}
    { split($2,a,":"); p=(a[3]a[4]=="AT"||a[3]a[4]=="TA"||a[3]a[4]=="CG"||a[3]a[4]=="GC") ? "palindromic" : "non-palindromic";
      m=$5; mb=(m<0.05)?"MAF<.05":(m<0.2)?"MAF.05-.2":"MAF>=.2"; r=d[$2]/N;
      for (k in K) delete K[k]; K[p]; K[mb]; K["chr"a[1]]; K["ALL"];
      for (k in K) { n[k]++; s[k]+=r; if (r>0.2) h[k]++ } }
    END{ for (k in n) printf "%-16s SNPs %8d  mean discordance %.4f  SNPs >20%% discordant %7d\n", k, n[k], s[k]/n[k], h[k] }' \
    cc_${2}_snp.txt cc_share_$2.frq | sort -k1,1V
  # per-woman discordance, by batch of the PAGE and SHARE copy
  awk 'NR>1{n[$3]++} END{for(i in n) print i, n[i]}' cc_$2.diff > cc_${2}_ind.txt
  awk -v V=$(wc -l < cc_share_$2.bim) 'FILENAME==ARGV[1]{b[$1" "$2]=$3; next} FILENAME==ARGV[2]{d[$1]=$2; next}
    { k=b[$2" P"]"/"b[$2" S"]; r=d[$2]/V; n[k]++; s[k]+=r; if(r>0.15) h[k]++ }
    END{ for(k in n) printf "copies %-18s women %4d  mean discordance %.4f  women >15%% %4d\n", k, n[k], s[k]/n[k], h[k] }' \
    dup_copy_batch.txt cc_${2}_ind.txt share_dup_rows.txt
}
conc "--extract prune.prune.in" pruned
conc "--chr 22" chr22
conc "--chr 2 --from-bp 53659707 --to-bp 55659707" sptbn1_region
plink --memory 60000 --bfile cc_share_sptbn1_region --bmerge cc_page_sptbn1_region --merge-mode 6 --snp 2:54659707:C:G --out cc_sptbn1 | grep -E 'overlapping|concordan' || true
echo "Done."
