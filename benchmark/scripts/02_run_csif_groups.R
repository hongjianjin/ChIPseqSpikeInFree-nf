# Step 2: ChIPseqSpikeInFree on the binned bigWig counts, one run per (cell system x mark), as in the paper.
# Uses the package's own functions (sourced from the package repository, no Bioconductor needed for these) in the order
# of the CLI's standard mode: ParseReadCounts -> CalculateSF -> PlotDistr -> BoxplotSF. Settings are the package defaults.
suppressMessages({library(parallel); library(readr)})
source("csif_pkg_ChIPseqSpikeInFree.R")   # = ChIPseqSpikeInFree/R/ChIPseqSpikeInFree.R, version 1.2.6
t1 <- read.csv("dataset1_table1_bw.csv", stringsAsFactors = FALSE)
t1$system <- sub("\\..*", "", t1$group); t1$cond <- sub(".*\\.", "", t1$group)
t1$run <- paste(t1$system, t1$mark, sep = "_")
dir.create("bench/csif", recursive = TRUE, showWarnings = FALSE)
run_group <- function(run) {
  d <- t1[t1$run == run, ]
  prefix <- file.path("bench/csif", run)
  first <- read_tsv(sprintf("bench/counts/%s_rawCounts.txt", d$gsm[1]), col_types = cols(bin = col_character(), .default = col_double()), progress = FALSE)
  vals <- list(first[[2]])
  for (g in d$gsm[-1]) {
    v <- read_tsv(sprintf("bench/counts/%s_rawCounts.txt", g), col_types = cols(bin = col_skip(), .default = col_double()), progress = FALSE)
    stopifnot(nrow(v) == nrow(first)); vals[[length(vals) + 1]] <- v[[1]]
  }
  df <- as.data.frame(vals, col.names = paste0(d$gsm, ".bw"), check.names = FALSE); rownames(df) <- first$bin
  meta <- data.frame(ID = paste0(d$gsm, ".bw"), ANTIBODY = d$mark, GROUP = d$cond)
  metaFile <- paste0(prefix, "_sample_meta.txt"); write.table(meta, metaFile, sep = "\t", quote = FALSE, row.names = FALSE)
  parsed <- ParseReadCounts(data = df, metaFile = metaFile, prefix = prefix, by = 0.05, binSize = 1000, ncores = 2)
  SF <- CalculateSF(parsed, metaFile = ReadMeta(metaFile), minFirstTurn = "auto", maxLastTurn = 0.99, cutoff_QC = 1.2)
  PlotDistr(parsed, SF, prefix, xlimMaxCPMW = NULL); BoxplotSF(SF, prefix)
  SF$gsm <- sub("\\.bw$", "", SF$ID); SF$run <- run; SF
}
res <- lapply(unique(t1$run), run_group)
out <- do.call(rbind, res); rownames(out) <- NULL
write.csv(out[, c("run", "gsm", "ID", "ANTIBODY", "GROUP", "QC", "TURNS", "SF")], "bench/dataset1_bigwig_csif_SF.csv", row.names = FALSE)
