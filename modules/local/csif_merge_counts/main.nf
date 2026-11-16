process CSIF_MERGE_COUNTS {
    tag "${meta.id}"
    label 'process_medium'

    container "${params.chipseqspikeinfree_container}"

    input:
    val  meta
    path counts, stageAs: 'counts/*'

    output:
    path "${prefix}_rawCounts.txt", emit: counts
    path "${prefix}_merge.log"    , emit: log
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat <<'END_R' > merge_counts.R
    args   <- commandArgs(trailingOnly = TRUE)
    prefix <- args[1]
    files  <- sort(args[-1])
    options(scipen = 99)

    # Samples dropped by CountRawReads (library size < 1e6) have no count column
    tabs <- lapply(files, function(f) {
        hdr <- strsplit(readLines(f, n = 1), "\\t", fixed = TRUE)[[1]]
        if (length(hdr) < 2 || !identical(hdr[1], "bin")) return(NULL)
        read.delim(f, check.names = FALSE, stringsAsFactors = FALSE,
                   colClasses = c("character", rep("integer", length(hdr) - 1)))
    })
    ok <- !vapply(tabs, is.null, logical(1))
    if (any(!ok)) {
        message("WARNING: dropped (library size < 1e6 reads or no counts): ",
                paste(basename(files[!ok]), collapse = ", "))
    }
    tabs <- tabs[ok]
    if (length(tabs) == 0) stop("No BAM has counts to merge")

    # Every BAM was counted on the same genome-wide bins, in the same order. Labels differ only
    # in chromosome notation (chr1 vs 1) for BAMs without a 'chr' prefix.
    norm <- function(x) sub("^M:", "MT:", sub("^chr", "", x))
    ref  <- norm(tabs[[1]][["bin"]])
    for (i in seq_along(tabs)) {
        if (nrow(tabs[[i]]) != length(ref) || !identical(norm(tabs[[i]][["bin"]]), ref)) {
            stop("Bins differ between per-BAM count files (check genome / bin size): ", files[ok][i])
        }
    }

    # As in CountRawReads: use the 'chr' labels when any BAM has them
    has_chr <- vapply(tabs, function(d) startsWith(d[["bin"]][1], "chr"), logical(1))
    bins    <- tabs[[if (any(has_chr)) which(has_chr)[1] else 1]][["bin"]]

    counts <- do.call(cbind, lapply(tabs, function(d) d[, -1, drop = FALSE]))
    out <- data.frame(bin = bins, counts, check.names = FALSE)
    write.table(out, paste0(prefix, "_rawCounts.txt"), sep = "\\t", quote = FALSE,
                row.names = FALSE, col.names = TRUE)
    message("Merged ", ncol(counts), " samples x ", nrow(out), " bins into ", prefix, "_rawCounts.txt")
    END_R

    Rscript merge_counts.R ${prefix} counts/* 2>&1 | tee ${prefix}_merge.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(paste(R.version\$major, R.version\$minor, sep = "."))' 2> /dev/null | tail -n 1)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_rawCounts.txt ${prefix}_merge.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: stub
    END_VERSIONS
    """
}
