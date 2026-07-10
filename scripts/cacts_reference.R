#!/usr/bin/env Rscript
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Harlan Barker
# Reference runner: score a representative matrix with the ORIGINAL CaCTS R (lawrenson-lab, GPL),
# for validating and benchmarking the pyCaCTS Python port. Sources their JSD.R from a LOCAL path
# (their GPL source is not redistributed in this repo). Loops every group as the query (their
# run_CaCTS_score scores one query at a time), records elapsed time, and writes the full score matrix.
# Usage: Rscript cacts_reference.R <rep_matrix.tsv> <out_scores.tsv> <path/to/JSD.R>
args <- commandArgs(trailingOnly = TRUE)
rep_path <- args[[1]]; out_path <- args[[2]]; jsd_r <- args[[3]]
suppressMessages(source(jsd_r))                                # provides run_CaCTS_score(), JSD(), KLD()
rep <- as.matrix(read.table(rep_path, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE))
dir.create("CaCTs_res", showWarnings = FALSE)                 # their fn writes a file per query
scores <- matrix(NA_real_, nrow(rep), ncol(rep), dimnames = list(rownames(rep), colnames(rep)))
t0 <- proc.time()[["elapsed"]]
invisible(capture.output(suppressMessages({
  for (q in colnames(rep)) {
    s <- run_CaCTS_score(rep, q)                              # df: Name, value(JSD), LogValue
    scores[as.character(s$Name), q] <- s$value
  }
})))
dt <- proc.time()[["elapsed"]] - t0
cat(sprintf("R_ELAPSED_SECONDS %.6f\n", dt))
write.table(data.frame(tf = rownames(scores), scores, check.names = FALSE),
            out_path, sep = "\t", quote = FALSE, row.names = FALSE)
unlink("CaCTs_res", recursive = TRUE)
