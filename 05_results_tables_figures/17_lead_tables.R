# 2026-11-25. Lead-SNP tables (Tables 1-2 and S5-S8 drafts): group-QC lead SNPs + Wald effect sizes (score-test p for significance).
# Effects per allele1 copy: beta in g/cm2, tau = log within-subject variance. Summary statistics only.
suppressMessages(library(data.table))
B <- file.path(Sys.getenv("WHI_DIR"), "10_Rebuild_v2/trajgwas"); S <- file.path(B, "summary")
sets <- c("black_hipbmd", "black_spinebmd", "hispanic_hipbmd", "hispanic_spinebmd", "ai_an_hipbmd")
out <- list()
for (s in sets) {
  ld <- fread(file.path(S, paste0("leads_", s, ".csv")), colClasses = list(character = c("snpid", "allele1", "allele2")))
  w <- list.files(file.path(B, s), pattern = "^chr[0-9]+\\.wald\\.txt$", full.names = TRUE)
  W <- rbindlist(lapply(w, fread, colClasses = list(character = c("snpid", "allele1", "allele2"))), fill = TRUE)
  W <- unique(W[snpid %in% ld$snpid, .(snpid, betaeffect, betastderr, wald_betap = betapval, taueffect, taustderr, wald_taup = taupval)], by = "snpid")
  x <- W[ld, on = "snpid"]
  x[, flag := fifelse(is.na(betaeffect), "no Wald", fifelse(!is.finite(betastderr) | !is.finite(taustderr) | taustderr > 5 | betastderr > 1, "check", ""))]
  out[[s]] <- x
}
TT <- rbindlist(out, fill = TRUE)
setcolorder(TT, c("set", "test", "chr", "pos", "snpid", "allele1", "allele2", "maf_grp", "mac", "hwe_grp", "betapval", "betaeffect", "betastderr", "wald_betap", "taupval", "taueffect", "taustderr", "wald_taup", "flag"))
fwrite(TT, file.path(S, "lead_snps_wald.csv"))
options(width = 250)
cat("Lead SNPs with Wald effects (score p = betapval/taupval; effects per allele1 copy). Wald found for", sum(TT$flag != "no Wald"), "of", nrow(TT), "; flagged 'check':", sum(TT$flag == "check"), "\n")
for (s in sets) { cat("\n==", s, "\n"); print(TT[set == s, .(test, snpid, maf = round(maf_grp, 3), p_beta = signif(betapval, 3), beta = signif(betaeffect, 3), se_b = signif(betastderr, 2),
                                                        p_tau = signif(taupval, 3), tau = signif(taueffect, 3), se_t = signif(taustderr, 2), flag)]) }
cat("\nFull table:", file.path(S, "lead_snps_wald.csv"), "\n")
