# 2026-11-25. FHS look-up summary: rebuilt-WHI Black lead SNPs + rs4233949 vs FHS women (same TrajGWAS model).
# WHI hip leads are looked up in FHS femoral neck and trochanter, spine leads in lumbar spine L2-4.
# TrajGWAS Wald effects are per copy of allele2 (checked 2026-11-25, 19_allele_check.R). FHS effects are first aligned to the
# WHI alleles, then both WHI and FHS effects are reported per copy of the WHI allele1 (sign flipped). Summary statistics only.
suppressMessages(library(data.table))
B <- file.path(Sys.getenv("WHI_DIR"), "10_Rebuild_v2"); TG <- file.path(B, "trajgwas"); S <- file.path(TG, "summary")
O <- file.path(Sys.getenv("FHS_DIR"), "10_Lookup_v2")
key <- function(chr, pos, a1, a2) paste(as.integer(chr), pos, pmin(a1, a2), pmax(a1, a2), sep = ":")
rd  <- function(f) if (file.exists(f)) fread(f, colClasses = list(character = c("snpid", "allele1", "allele2"))) else NULL
W <- list()
for (s in c("black_hipbmd", "black_spinebmd")) {
  ld <- fread(file.path(S, paste0("leads_", s, ".csv")), colClasses = list(character = c("snpid", "allele1", "allele2")))
  ld <- ld[, .(set = s, lead = ifelse(grepl("beta", test), "mean", "variability"), chr, pos, snpid, allele1, allele2, maf_whi = maf_grp, whi_score_betap = betapval, whi_score_taup = taupval)]
  sp <- rd(Sys.glob(file.path(TG, s, "*sptbn1*wald*.txt"))[1])
  if (!is.null(sp)) ld <- rbind(ld, sp[1, .(set = s, lead = "rs4233949", chr, pos, snpid, allele1, allele2, maf_whi = if ("maf" %in% names(sp)) maf else NA_real_, whi_score_betap = NA_real_, whi_score_taup = NA_real_)])
  wd <- rbindlist(lapply(Sys.glob(file.path(TG, s, "chr*.wald.txt")), rd), fill = TRUE)
  if (!is.null(sp)) wd <- rbind(wd, sp, fill = TRUE)
  if (nrow(wd)) ld <- merge(ld, unique(wd[, .(snpid, whi_beta = betaeffect, whi_beta_se = betastderr, whi_tau = taueffect, whi_tau_se = taustderr)], by = "snpid"), by = "snpid", all.x = TRUE)
  W[[s]] <- ld
}
W <- rbindlist(W, fill = TRUE); W[, k := key(chr, pos, allele1, allele2)]
match_o <- list(black_hipbmd = c("fnbmd", "trbmd"), black_spinebmd = "lsbmd")
R <- list()
for (o in c("fnbmd", "trbmd", "lsbmd")) {
  sc <- rd(file.path(O, paste0("fhs_", o, ".score.txt"))); wa <- rd(file.path(O, paste0("fhs_", o, ".wald.txt")))
  if (is.null(sc) || is.null(wa)) { cat("no FHS results for", o, "\n"); next }
  fq <- fread(file.path(O, paste0("lookup_", o, ".frq")))[, .(snpid = SNP, maf_fhs = MAF, mac_fhs = round(MAF * NCHROBS))]
  f <- merge(sc[, .(snpid, chr, pos, f_a1 = allele1, f_a2 = allele2, fhs_score_betap = betapval, fhs_score_taup = taupval)],
             wa[, .(snpid, fhs_beta = betaeffect, fhs_beta_se = betastderr, fhs_tau = taueffect, fhs_tau_se = taustderr)], by = "snpid")
  f <- merge(f, fq, by = "snpid"); f[, k := key(chr, pos, f_a1, f_a2)]
  n_women <- length(readLines(file.path(O, paste0("lookup_", o, ".fam"))))
  w <- W[vapply(set, function(s) o %in% match_o[[s]], TRUE)]
  x <- merge(w, f[, !c("chr", "pos", "snpid")], by = "k", all.x = TRUE)
  x[, flip := !is.na(f_a1) & f_a1 != allele1]
  x[flip == TRUE, `:=`(fhs_beta = -fhs_beta, fhs_tau = -fhs_tau)]
  x[, `:=`(fhs_outcome = o, fhs_women = n_women, fhs_testable = !is.na(maf_fhs) & maf_fhs >= 0.01 & mac_fhs >= 20)]
  R[[o]] <- x
}
R <- rbindlist(R, fill = TRUE)
for (v in c("whi_beta", "whi_tau", "fhs_beta", "fhs_tau")) set(R, j = v, value = -R[[v]])   # per copy of WHI allele1
R[, same_dir := fifelse(lead == "variability", sign(fhs_tau) == sign(whi_tau), sign(fhs_beta) == sign(whi_beta))]
R[, fhs_p := fifelse(lead == "variability", fhs_score_taup, fhs_score_betap)]
out <- R[, .(set, lead, fhs_outcome, fhs_women, snpid, allele1, allele2, maf_whi, maf_fhs, mac_fhs, fhs_testable,
             whi_score_betap, whi_beta, whi_beta_se, fhs_score_betap, fhs_beta, fhs_beta_se,
             whi_score_taup, whi_tau, whi_tau_se, fhs_score_taup, fhs_tau, fhs_tau_se, same_dir)]
setorder(out, set, fhs_outcome, lead)
fwrite(out, file.path(O, "fhs_lookup_results.csv")); fwrite(out, file.path(S, "fhs_lookup_results.csv"))
options(width = 250)
cat("\nFHS women (unrelated, with outcome):", paste(unique(paste(R$fhs_outcome, R$fhs_women)), collapse = ", "), "\n")
cat("\nrs4233949 (effects per copy of WHI allele1):\n")
print(out[lead == "rs4233949", .(set, fhs_outcome, allele1, maf_whi, maf_fhs, whi_beta, whi_beta_se, fhs_beta, fhs_beta_se, fhs_score_betap, fhs_tau, fhs_score_taup)], digits = 3)
cat("\nLead SNPs by discovery set x FHS outcome (tested = MAF >= .01 and MAC >= 20 in FHS):\n")
out[, fhs_p := fifelse(lead == "variability", fhs_score_taup, fhs_score_betap)]
print(out[lead != "rs4233949", .(n_leads = .N, tested = sum(fhs_testable), same_direction = sum(same_dir & fhs_testable, na.rm = TRUE),
         p_lt_05 = sum(fhs_testable & fhs_p < 0.05, na.rm = TRUE), bonferroni = signif(0.05 / max(1, sum(fhs_testable)), 2),
         min_p = suppressWarnings(min(fhs_p[fhs_testable], na.rm = TRUE))), by = .(set, lead, fhs_outcome)], digits = 3)
cat("\nLead SNPs with FHS p < 0.05 (test matching the discovery test):\n")
print(out[lead != "rs4233949" & fhs_testable & fhs_p < 0.05,
          .(set, lead, fhs_outcome, snpid, maf_whi, maf_fhs, whi_beta, fhs_beta, whi_tau, fhs_tau, fhs_score_betap, fhs_score_taup, same_dir)], digits = 3)
cat("\nFull table:", file.path(S, "fhs_lookup_results.csv"), "\n")
