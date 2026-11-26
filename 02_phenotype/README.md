# Phenotype files

`build_pheno_files.R` builds the long-format phenotype files read by TrajGWAS directly from the WHI dbGaP tables
(one row per woman per DXA visit):

| Column | Source |
|---|---|
| `time` | DXA visit year (`HIPVY` / `SPNVY`) |
| `hipbmd` | `HIPBMDCOR`, total hip BMD, scanner-corrected (`bmdhip_rel2`) |
| `spinebmd` | `SPNBMDCOR`, total spine BMD, scanner-corrected (`bmdspine_rel2`) |
| `AGE`, `RACE` | baseline age and race/ethnicity (`f2_rel1`) |
| `HEIGHTX`, `WEIGHTX` | baseline height and weight, visit year 0 (`f80_rel1`) |

Genotype IDs are linked to dbGaP subjects through the WHI Sample file (`SAMPLE_ID`) and Subject file
(`SUBJID`; `SOURCE_SUBJECT_ID` for GARNET). Usage: `source config.sh; Rscript 02_phenotype/build_pheno_files.R [fam] [out]`.

`check_bmd_source.R` confirms that the BMD values equal the scanner-corrected dbGaP variables.
