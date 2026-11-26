#!/bin/bash
#SBATCH --job-name=page_reimp_dl
#SBATCH --output=page_reimp_dl_%A_%a.out
#SBATCH --time=12:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=4
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --mail-type=FAIL
#SBATCH --array=1-9
# 2026-11-24. Download + convert the 9 PAGE chromosomes re-imputed on Michigan (after removing
# the 9 low-call-rate women): c1 chr6, 9, 17, 19 and c2 chr6, 14, 20, 21, 22.
# Before submitting:
#   ~/mis_dl_page_c1.txt, ~/mis_dl_page_c2.txt : Michigan download commands for each job
#   ~/.mis_pw_page_c1,    ~/.mis_pw_page_c2    : the result password (chmod 600)
# Variant IDs are rebuilt as chr:pos:REF:ALT so they match the rest of PAGE_cX_R2.
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # site paths (see config.sh)
TASKS=(c1:6 c1:9 c1:17 c1:19 c2:6 c2:14 c2:20 c2:21 c2:22)
T=${TASKS[$((SLURM_ARRAY_TASK_ID-1))]}; CG=${T%%:*}; i=${T##*:}
OUT=${WHI_DBGAP_DIR}/02_PAGE/05_Rebuild/reimpute/$CG/results
mkdir -p $OUT; cd $OUT
echo "PAGE $CG chr$i"

# 1. download (skip the "echo Downloading" lines in Michigan's list)
CMD=$(grep -E "chr_${i}\.zip" ~/mis_dl_page_${CG}.txt | grep -vE "^[[:space:]]*echo" | head -1)
[ -n "$CMD" ] || { echo "No download line for chr_${i}.zip in ~/mis_dl_page_${CG}.txt"; exit 1; }
eval "$CMD"
[ -s chr_${i}.zip ] || { echo "Download failed"; exit 1; }

# 2. extract dose + info
ok=0
while IFS= read -r PW; do
  if 7z x -y -p"$PW" chr_${i}.zip "chr${i}.dose.vcf.gz" "chr${i}.info.gz" >/dev/null 2>&1; then ok=1; break; fi
done < ~/.mis_pw_page_${CG}
[ $ok -eq 1 ] || { echo "Could not unzip with the password"; exit 1; }
mv chr_${i}.zip chr_${i}.zip.done   # kept until the rebuild is checked; delete afterwards

# 3. R2>=0.3 list as chr:pos:REF:ALT (Michigan 2 info is VCF-style, R2= in INFO)
zcat chr${i}.info.gz | awk -v OFS=: '/^#/{next}
  $8 ~ /R2=/ { c=$1; sub(/^chr/,"",c); if (match($8,/R2=[0-9.eE-]+/)) { r=substr($8,RSTART+3,RLENGTH-3); if (r+0>=0.3) print c,$2,$4,$5 } }' > r2pass_chr${i}.txt

# 4. convert with REF/ALT kept (A2=REF, A1=ALT), rename IDs, then keep R2>=0.3
plink --memory 40000 --vcf chr${i}.dose.vcf.gz --double-id --keep-allele-order --make-bed --out tmp_chr${i} > /dev/null
awk -v OFS='\t' '{c=$1; sub(/^chr/,"",c); $1=c; $2=c":"$4":"$6":"$5; print}' tmp_chr${i}.bim > tmp_chr${i}.bim.new
mv tmp_chr${i}.bim.new tmp_chr${i}.bim
plink --memory 40000 --bfile tmp_chr${i} --keep-allele-order --extract r2pass_chr${i}.txt --make-bed --out chr${i} > /dev/null
echo "chr$i: imputed $(zcat chr${i}.info.gz | grep -vc '^#'), R2>=0.3 $(wc -l < r2pass_chr${i}.txt), kept $(wc -l < chr${i}.bim) variants, $(wc -l < chr${i}.fam) people; first ID $(head -1 chr${i}.bim | cut -f2)"
[ "$i" = 6 ] && [ "$CG" = c1 ] && echo "old PAGE_c1 chr6 R2>=0.3 was 654,100" || true
