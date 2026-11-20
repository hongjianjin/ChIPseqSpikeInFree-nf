# Step 3: how stable is CalculateSF? (a) recompute SF from the parsedMatrix files as CSIF_SCALE does (matrix read from file);
# (b) recompute with tiny relative noise on the cumulative values (1e-14, 1e-10; 50 replicates); (c) same with minFirstTurn fixed at 0.5.
suppressMessages(library(parallel)); source("csif_pkg_ChIPseqSpikeInFree.R")
rd <- function(f) read.table(f, sep = "\t", header = TRUE, quote = "", check.names = FALSE, row.names = NULL)
runs <- c("DIPG_H3K27Ac", "DIPG_H3K27me3", "DIPG_H3K4me3", "H-NSC_H3K27me3")
set.seed(42); file_path <- list(); jit <- list()
for (run in runs) {
  M <- rd(sprintf("bench/csif/%s_parsedMatrix.txt", run)); meta <- ReadMeta(sprintf("bench/csif/%s_sample_meta.txt", run))
  s0 <- CalculateSF(M, meta)
  file_path[[run]] <- data.frame(run = run, gsm = sub("\\.bw$", "", s0$ID), sf_file = s0$SF, QC_file = s0$QC, TURNS_file = s0$TURNS)
  for (eps in c(1e-14, 1e-10)) for (rep in 1:50) {
    Mp <- M; Mp[, -1] <- as.matrix(M[, -1]) * (1 + eps * runif(nrow(M) * (ncol(M) - 1), -1, 1))
    s <- CalculateSF(Mp, meta)
    jit[[length(jit) + 1]] <- data.frame(run = run, gsm = sub("\\.bw$", "", s$ID), eps = eps, rep = rep, sf = s$SF, qc = s$QC)
  }
}
write.csv(do.call(rbind, file_path), "bench/dataset1_bigwig_csif_SF_from_file.csv", row.names = FALSE)
write.csv(do.call(rbind, jit), "bench/dataset1_bigwig_csif_SF_jitter_long.csv", row.names = FALSE)
# (c) exploratory: first turning point fixed (code path without the automatic density step)
set.seed(5); out <- list()
for (run in runs) {
  M <- rd(sprintf("bench/csif/%s_parsedMatrix.txt", run)); meta <- ReadMeta(sprintf("bench/csif/%s_sample_meta.txt", run))
  for (mf in c("auto", "0.5")) {
    fm <- if (mf == "auto") "auto" else as.numeric(mf)
    base <- CalculateSF(M, meta, minFirstTurn = fm)
    reps <- sapply(1:50, function(i) { Mp <- M; Mp[, -1] <- as.matrix(M[, -1]) * (1 + 1e-14 * runif(nrow(M) * (ncol(M) - 1), -1, 1)); CalculateSF(Mp, meta, minFirstTurn = fm)$SF })
    out[[length(out) + 1]] <- data.frame(run = run, minFirstTurn = mf, samples = nrow(meta), unstable = sum(apply(reps, 1, function(v) length(unique(v)) > 1)),
                                         gsm = I(list(sub("\\.bw$", "", base$ID))), sf = I(list(base$SF)))
  }
}
res <- do.call(rbind, out); print(res[, c("run", "minFirstTurn", "samples", "unstable")])
fixed <- do.call(rbind, lapply(seq_len(nrow(res)), function(i) if (res$minFirstTurn[i] == "0.5") data.frame(run = res$run[i], gsm = res$gsm[[i]], `sf_fixed0.5` = res$sf[[i]], check.names = FALSE)))
write.csv(fixed, "bench/dataset1_bigwig_csif_SF_minFirstTurn0.5.csv", row.names = FALSE)
