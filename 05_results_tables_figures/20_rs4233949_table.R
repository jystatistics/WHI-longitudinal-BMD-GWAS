# 2026-11-25. rs4233949 (SPTBN1) in the rebuilt WHI (5 sets) and FHS women: score p (SPA) + Wald effects per copy of allele1.
# TrajGWAS Wald effects are per allele2, so signs are flipped. Summary statistics only.
suppressMessages(library(data.table))
B <- file.path(Sys.getenv("WHI_DIR"), "10_Rebuild_v2/trajgwas"); O <- file.path(Sys.getenv("FHS_DIR"), "10_Lookup_v2")
cc <- list(character = c("snpid", "allele1", "allele2"))
R <- list()
for (s in c("black_hipbmd", "black_spinebmd", "hispanic_hipbmd", "hispanic_spinebmd", "ai_an_hipbmd")) {
  sc <- fread(file.path(B, s, "chr2.sptbn1_rs4233949.pval.txt"), colClasses = cc)[1]
  wa <- fread(file.path(B, s, "chr2.sptbn1_rs4233949.wald.txt"), colClasses = cc)[1]
  frq <- file.path(B, "summary", paste0("grp_", s, ".frq"))
  R[[s]] <- data.table(study = "WHI", set = s, allele1 = sc$allele1, allele2 = sc$allele2, maf_all = sc$maf,
                       beta = -wa$betaeffect, beta_se = wa$betastderr, p_beta = sc$betapval, tau = -wa$taueffect, tau_se = wa$taustderr, p_tau = sc$taupval)
}
for (o in c("fnbmd", "lsbmd")) {
  sc <- fread(file.path(O, paste0("fhs_", o, ".score.txt")), colClasses = cc)[pos == 54659707]
  wa <- fread(file.path(O, paste0("fhs_", o, ".wald.txt")), colClasses = cc)[pos == 54659707]
  fq <- fread(file.path(O, paste0("lookup_", o, ".frq")))[grepl("^2:54659707", SNP)]
  sg <- if (sc$allele1 == "C") -1 else 1          # express per C allele (WHI allele1)
  R[[o]] <- data.table(study = "FHS", set = o, allele1 = "C", allele2 = "G", maf_all = fq$MAF,
                       beta = sg * wa$betaeffect, beta_se = wa$betastderr, p_beta = sc$betapval, tau = sg * wa$taueffect, tau_se = wa$taustderr, p_tau = sc$taupval)
}
R <- rbindlist(R)
# group-specific C-allele frequency for WHI sets (SPTBN1 is not in WHI_v5_r08 at R2 >= 0.8; use the look-up file, FID = IID)
L <- file.path(Sys.getenv("WHI_DIR"), "10_Rebuild_v2/lookup/sptbn1_rs4233949")
for (s in R[study == "WHI", set]) {
  k <- fread(file.path(B, "summary", paste0("keep_", s, ".txt")), header = FALSE, colClasses = "character")
  tk <- tempfile(); writeLines(paste(k$V2, k$V2), tk); to <- tempfile()
  system(sprintf("plink --bfile %s --keep %s --freq --out %s > /dev/null", L, tk, to))
  f <- tryCatch(fread(paste0(to, ".frq")), error = function(e) NULL)
  if (!is.null(f) && nrow(f)) R[set == s, `:=`(eaf = if (f$A1[1] == allele1) f$MAF[1] else 1 - f$MAF[1], n_geno = f$NCHROBS[1] / 2)]
}
R[study == "FHS", eaf := maf_all]
fwrite(R, file.path(B, "summary", "rs4233949_whi_fhs.csv"))
options(width = 200); print(R[, .(study, set, EA = allele1, eaf = round(eaf, 3), beta = signif(beta, 3), beta_se = signif(beta_se, 2), p_beta = signif(p_beta, 3),
                                  tau = signif(tau, 3), tau_se = signif(tau_se, 2), p_tau = signif(p_tau, 3))])
