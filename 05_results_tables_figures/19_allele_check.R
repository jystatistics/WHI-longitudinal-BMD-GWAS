# 2026-11-25. Which allele does the TrajGWAS Wald effect count? Compare the Wald beta sign with a plain regression of
# baseline-visit hip BMD on the PLINK A1 count for the top Black hip mean-BMD SNPs. Prints aggregate numbers only.
suppressMessages(library(data.table))
B <- file.path(Sys.getenv("WHI_DIR"), "10_Rebuild_v2"); S <- file.path(B, "trajgwas", "summary")
ld <- fread(file.path(S, "lead_snps_wald.csv"), colClasses = list(character = c("snpid", "allele1", "allele2")))[set == "black_hipbmd" & grepl("beta", test)][1:3]
writeLines(ld$snpid, file.path(S, "allele_check_snps.txt"))
system(sprintf("plink --bfile %s/WHI_v5_r08 --keep %s/keep_black_hipbmd.txt --extract %s/allele_check_snps.txt --recode A --out %s/allele_check > /dev/null", B, S, S, S))
g <- fread(file.path(S, "allele_check.raw")); unlink(file.path(S, "allele_check.raw"))
al <- fread(file.path(B, "alias_all.txt"), header = FALSE, colClasses = "character"); amap <- setNames(al$V2, al$V1)
ph <- fread(file.path(Sys.getenv("WHI_DIR"), "01_Black/pheno_black_hipbmd_data.csv"))
ph[, IID := as.character(IID)]; j <- !is.na(amap[ph$IID]); ph$IID[j] <- amap[ph$IID[j]]
ph <- ph[is.finite(hipbmd)][order(IID, time)][, .SD[1], by = IID]          # first visit per woman
g[, IID := as.character(IID)]; d <- merge(ph[, .(IID, hipbmd)], g, by = "IID")
bim <- fread(file.path(B, "WHI_v5_r08.bim"), header = FALSE)[V2 %in% ld$snpid, .(snpid = V2, plinkA1 = V5, plinkA2 = V6)]
for (i in seq_len(nrow(ld))) {
  col <- grep(paste0("^", ld$snpid[i], "_"), names(d), value = TRUE, fixed = FALSE)
  cnt <- sub(".*_", "", col); x <- d[[col]]
  fit <- coef(summary(lm(d$hipbmd ~ x)))
  cat(sprintf("%s | TrajGWAS allele1=%s allele2=%s, Wald beta %+.4f | plink A1=%s, counted allele in .raw=%s, crude slope per copy %+.4f (n=%d); mean BMD by count: %s\n",
      ld$snpid[i], ld$allele1[i], ld$allele2[i], ld$betaeffect[i], bim[snpid == ld$snpid[i], plinkA1], cnt, fit[2, 1], sum(!is.na(x)),
      paste(tapply(d$hipbmd, x, function(v) sprintf("%.3f (n=%d)", mean(v), length(v))), collapse = ", ")))
}
cat("If the crude slope has the same sign as the Wald beta, the Wald effect is per copy of the counted .raw allele; if opposite, per the other allele.\n")
