# 2026-11-25. Pick SNPs for the TrajGWAS Wald test (effect sizes): per set, the top 10 by the mean (beta) test,
# top 10 by the variability (tau) test, top 10 by the joint test, and the best beta/tau/joint SNP in the SPTBN1 region.
suppressMessages(library(data.table))
B <- file.path(Sys.getenv("WHI_DIR"), "10_Rebuild_v2/trajgwas")
sets <- c("black_hipbmd", "black_spinebmd", "hispanic_hipbmd", "hispanic_spinebmd", "ai_an_hipbmd")
for (s in sets) {
  d <- rbindlist(lapply(file.path(B, s, sprintf("chr%d.pval.txt", 1:22)), fread))
  for (v in c("maf", "hwepval", "betapval", "taupval", "jointpval")) set(d, j = v, value = suppressWarnings(as.numeric(d[[v]])))
  d <- d[maf >= 0.01 & hwepval >= 1e-6 & is.finite(jointpval)]
  r <- d[chr == 2 & pos >= 54400000 & pos <= 55200000]
  ids <- unique(c(head(d[order(betapval)], 10)$snpid, head(d[order(taupval)], 10)$snpid, head(d[order(jointpval)], 10)$snpid,
                  r[which.min(betapval)]$snpid, r[which.min(taupval)]$snpid, r[which.min(jointpval)]$snpid))
  writeLines(ids, file.path(B, "summary", paste0("wald_snps_", s, ".txt")))
  cat(sprintf("%s: %d SNPs on chr %s\n", s, length(ids), paste(sort(unique(as.integer(sub(":.*", "", ids)))), collapse = ",")))
  rm(d); gc()
}
