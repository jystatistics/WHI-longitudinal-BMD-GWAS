# 2026-11-25. FHS look-up, step 1: women with the outcome and covariates who are in the FHS female genotype file.
# Args: phenotype csv, outcome column, fam file, output keep list. Prints counts only.
a <- commandArgs(TRUE); suppressMessages(library(data.table))
ph <- fread(a[1]); y <- a[2]; ph[, IID := as.character(IID)]
for (v in c(y, "time", "AGE20", "BMI20")) set(ph, j = v, value = suppressWarnings(as.numeric(ph[[v]])))
cat("phenotype rows", nrow(ph), "| sex values:", paste(names(table(ph$sex)), table(ph$sex), sep = "=", collapse = " "), "\n")
ph <- ph[is.finite(get(y)) & is.finite(time) & is.finite(AGE20) & is.finite(BMI20)]
fl <- strsplit(trimws(readLines(a[3])), "[[:space:]]+"); fam <- data.table(FID = vapply(fl, `[`, "", 1), IID = vapply(fl, `[`, "", 2))
cat("fam: ", nrow(fam), "people; FID == IID for", sum(fam$FID == fam$IID), "| unique IID", uniqueN(fam$IID), "\n")
ph <- ph[IID %in% fam$IID]
nv <- ph[, .N, by = IID]$N
cat("complete rows in genotype file:", nrow(ph), "visits,", uniqueN(ph$IID), "women; visits per woman:",
    paste(names(table(nv)), table(nv), sep = "=", collapse = ", "), "\n")
cat("time: range", paste(round(range(ph$time), 2), collapse = " to "), "| most common values:",
    paste(head(names(sort(table(round(ph$time, 1)), decreasing = TRUE)), 10), collapse = " "), "\n")
cat("AGE20 constant within woman:", ph[, uniqueN(AGE20), by = IID][, all(V1 == 1)], "| range", paste(round(range(ph$AGE20), 1), collapse = "-"),
    "| BMI20 constant:", ph[, uniqueN(BMI20), by = IID][, all(V1 == 1)], "| outcome mean", round(mean(ph[[y]]), 3), "sd", round(sd(ph[[y]]), 3), "\n")
k <- fam[IID %in% ph$IID]
writeLines(paste(k$FID, k$IID), a[4]); cat("keep list:", nrow(k), "women ->", a[4], "\n")
