# 2026-11-25. Combine TrajGWAS Wald results (effect sizes) with the score-test p-values for the same SNPs.
suppressMessages(library(data.table))
B <- file.path(Sys.getenv("WHI_DIR"), "10_Rebuild_v2/trajgwas")
sets <- c("black_hipbmd", "black_spinebmd", "hispanic_hipbmd", "hispanic_spinebmd", "ai_an_hipbmd")
out <- list()
for (s in sets) {
  w <- list.files(file.path(B, s), pattern = "^chr[0-9]+\\.wald\\.txt$", full.names = TRUE)
  if (!length(w)) next
  W <- rbindlist(lapply(w, fread))
  chrs <- unique(W$chr)
  S <- rbindlist(lapply(file.path(B, s, sprintf("chr%d.pval.txt", chrs)), fread))[snpid %in% W$snpid, .(snpid, score_betapval = betapval, score_taupval = taupval)]
  W <- S[W, on = "snpid"][, set := s]
  # flag fits that look unreliable: non-finite or very large SEs
  W[, flag := fifelse(!is.finite(betastderr) | !is.finite(taustderr) | taustderr > 5 | betastderr > 1, "check", "")]
  out[[s]] <- W
}
TT <- rbindlist(out, fill = TRUE)
setcolorder(TT, c("set", "chr", "pos", "snpid", "allele1", "allele2", "maf", "betaeffect", "betastderr", "betapval", "score_betapval", "taueffect", "taustderr", "taupval", "score_taupval", "flag"))
fwrite(TT, file.path(B, "summary", "trajgwas_wald_effects.csv"))
options(width = 250)
cat("Wald effect sizes (beta = BMD g/cm2 per allele1 copy; tau = log within-subject variance per copy). Rows flagged 'check' have non-finite/very large SEs.\n")
print(TT[order(set, pmin(betapval, taupval))][, head(.SD, 6), by = set][, .(set, snpid, maf, betaeffect, betastderr, betapval, score_betapval, taueffect, taustderr, taupval, score_taupval, flag)], digits = 3)
cat(sprintf("\n%d rows; %d flagged 'check'. Full table: %s\n", nrow(TT), sum(TT$flag == "check"), file.path(B, "summary", "trajgwas_wald_effects.csv")))
