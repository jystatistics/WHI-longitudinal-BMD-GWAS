# Longitudinal GWAS of bone mineral density in the Women's Health Initiative

Analysis code for:

> Jung J, Wu Q. *Longitudinal Genome-Wide Association Study of Bone Mineral Density in Black, Hispanic, and
> American Indian/Alaska Native Postmenopausal Women from the Women's Health Initiative.* PLOS Aging and Health
> (under revision, manuscript PAGE-D-26-00166).

This repository contains code only. It holds **no individual-level data, participant IDs, phenotype files,
or access credentials**. The data are available to approved users from dbGaP (see *Data access*).

## What the pipeline does

1. **Genotype rebuild (`01_genotype_rebuild/`).** Re-imputes and merges five WHI genotype studies
   (GARNET, SHARE, PAGE, WHIMS, WHISP) on GRCh37 with one pipeline. Variants are renamed `chr:pos:REF:ALT`,
   kept only if imputation R² ≥ 0.8 in every batch imputed for that chromosome, and women genotyped by more
   than one study are kept once. Output: `WHI_v5_r08` (2,363 unique women) and genetic principal components
   (`WHI_PCA_v5_r08.eigenvec`).
2. **Phenotypes (`02_phenotype/`).** Builds one long-format file per group and site from the WHI dbGaP tables:
   scanner-corrected total hip and spine BMD (`HIPBMDCOR`, `SPNBMDCOR`) at each DXA visit, with baseline age,
   height, weight and race.
3. **TrajGWAS (`03_trajgwas/`).** One pre-specified model for every group and skeletal site:
   * mean model: BMD ~ 1 + time + age + height + weight + WHI sub-study + PC1–PC10
   * within-subject variance model: ~ 1 + age + BMI
   * random effects: intercept and slope for time (years since the baseline DXA visit)
   * score tests genome-wide; Wald tests for effect sizes at lead variants and rs4233949
4. **FHS look-up (`04_fhs_lookup/`).** The same model in unrelated Framingham Heart Study women, for the
   femoral neck and lumbar spine, at the WHI lead variants and rs4233949.
5. **Results (`05_results_tables_figures/`).** Manhattan and QQ plots, lead-variant tables, gene annotation,
   allele checks, and the rs4233949 table.

## Run order

All jobs were run on a SLURM cluster. Set the paths in `config.sh`, replace `YOUR_ACCOUNT` in the
`#SBATCH` headers, and submit each job from the directory that contains `config.sh`.

| Step | Script | Purpose |
|---|---|---|
| 1 | `01_genotype_rebuild/01_SHARE_prep_for_TOPMed.sh` | SHARE: lift hg18 → hg19, align to HRC, write VCFs for imputation |
| 2 | `01_genotype_rebuild/01b_SHARE_rewrite_vcf_ref.sh` | Set REF/ALT to the HRC reference allele before upload |
| – | Michigan Imputation Server | HRC r1.1 (GRCh37), Eagle v2.4 phasing, population "mixed", rsq filter off |
| 3 | `02_PAGE_convert_imputed.sh` | PAGE imputed dosages → PLINK |
| 4 | `03_SHARE_download_convert.sh`, `03b_SHARE_merge_chr.sh`, `03c_SHARE_fixid_merge.sh` | SHARE imputed results → PLINK, `chr:pos:REF:ALT` IDs |
| 5 | `05_PAGE_reimpute_prep.sh`, `05b_PAGE_reimpute_refalign.sh` | PAGE chromosomes whose chunks the server excluded: drop low-call-rate samples, REF-align, re-upload |
| 6 | `08_PAGE_reimputed_download.sh`, `08b_PAGE_swap_chr.sh` | Bring the re-imputed PAGE chromosomes back in |
| 7 | `04_merge_all_cohorts.sh` | Merge all batches (WHISP genotypes set to missing where WHISP fails the R² rule) |
| 8 | `04b_dedup_page_share.sh`, `04c_dup_concordance.sh` | Remove women present in both PAGE and SHARE; genotype concordance of the duplicate copies |
| 9 | `04d_r2_strict_rebuild.sh` | R² ≥ 0.8 rule, LD-pruned PCA |
| 10 | `04e_dedup_cross_id.sh`, `04f_dedup_z2.sh` | Remove women genotyped under different IDs (IBS/IBD), final `WHI_v5_r08` + PCs |
| 11 | `02_phenotype/build_pheno_files.R` | Phenotype files from the dbGaP tables (`check_bmd_source.R`: source check) |
| 12 | `03_trajgwas/06_trajgwas_setup.sh` | Split genotypes by chromosome, sub-study labels |
| 13 | `03_trajgwas/08_submit_trajgwas.sh full` | Genome-wide TrajGWAS (`07_trajgwas_run.sh` → `07_trajgwas.jl`), then `09_summary.sh` |
| 14 | `09b_trajgwas_calibration.R`, `09c_lambda_dirs.R` | Genomic inflation (λGC) |
| 15 | `08_submit_trajgwas.sh lookup` | rs4233949 look-up |
| 16 | `08_submit_trajgwas.sh leads` | Group-specific QC (`13_group_qc.R`), lead loci, Wald tests |
| 17 | `04_fhs_lookup/16a_fhs_prep.sh` → `16c_fhs_run.sh` → `16e_fhs_summary.R` | FHS look-up |
| 18 | `08_submit_trajgwas.sh figures`; `05_results_tables_figures/*.R` | Figures and tables |

`10_wald_select.R` and `10b_wald_table.R` are an earlier Wald step (top 10 variants per test) kept for
completeness; the published lead tables come from steps 16 and 18.

The scripts in `03_trajgwas/`, `04_fhs_lookup/` and `05_results_tables_figures/` were run from
`$WHI_DIR/10_Rebuild_v2`; copy them there (with `config.sh`) before submitting. `07_trajgwas.jl` also has
diagnostic modes (`DIAG=1`, `DIAG=7`) that we used to choose the optimizer settings; they do not change
the model.

## Software

* PLINK 1.9; bcftools 1.17; UCSC liftOver (hg18ToHg19 chain); 7-Zip
* R 4.4.0 with `data.table` and `ragg`
* Julia 1.6.5 with `TrajGWAS`, `WiSER`, `SnpArrays`, `CSV`, `DataFrames`
* Michigan Imputation Server 2, HRC r1.1 reference panel (GRCh37)
* `HRC.r1-1.GRCh37.wgs.mac5.sites.tab` (HRC site list) in `$TOOLS_DIR`

## Data access

Individual-level data cannot be redistributed under the dbGaP data use agreements. Approved users can
request them from dbGaP:

* WHI: phs000281, phs000675, phs000227, phs000315, phs000386
* Framingham Heart Study: phs000007

Michigan Imputation Server download links and result passwords are read from files in the user's home
directory (`~/mis_dl_*.txt`, `~/.mis_pw_*`, `chmod 600`) and are never stored in this repository.

## License

MIT (see `LICENSE`). If you use this code, please cite the article and the Zenodo archive (`CITATION.cff`).
