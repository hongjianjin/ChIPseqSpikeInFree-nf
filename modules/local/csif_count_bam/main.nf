process CSIF_COUNT_BAM {
    tag "${meta.id}"
    label 'process_low'
    label 'process_long'

    container "${params.chipseqspikeinfree_container}"

    input:
    tuple val(meta), path(bam, stageAs: 'bam/*')
    path  chrom_sizes, stageAs: 'chrom_sizes/*'
    val   genome
    val   bin_size
    val   paired_end

    output:
    tuple val(meta), path("${meta.id}_rawCounts.txt"), emit: counts
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def chrom      = chrom_sizes ? chrom_sizes : genome
    def single_end = paired_end ? 'FALSE' : 'TRUE'
    """
    # ChIPseqSpikeInFree names the count column after the BAM file: link it as <sample>.bam
    ln -s "${bam}" "${meta.id}.bam"

    cat <<'END_R' > count_bam.R
    # Counting as in ChIPseqSpikeInFree::CountRawReads(), for one BAM. CountRawReads() itself
    # cannot take a single BAM: it drops the bin column with datOut[, -1], gets a vector, and
    # fails in colSums(). The steps below are its own, applied to one BAM.
    suppressMessages({
        library(ChIPseqSpikeInFree)
        library(Rsamtools)
        library(GenomicAlignments)
    })
    args <- commandArgs(trailingOnly = TRUE)
    bam <- args[1]; chrom <- args[2]; prefix <- args[3]
    single_end <- as.logical(args[4]); bin_size <- as.numeric(args[5])
    min_reads <- 1e6   # CountRawReads drops samples with total counts not above 10^6
    out <- paste0(prefix, "_rawCounts.txt")

    # Bins named with or without 'chr', matching the BAM header (CountRawReads: ChrCheck)
    targets  <- names(scanBamHeader(bam)[[1]][["targets"]])
    with_chr <- any(grepl("^chr", targets))
    bins <- GenerateBins(chromFile = chrom, binSize = bin_size, overlap = 0, withChr = with_chr)
    coords <- GRanges(
        seqnames = bins[, 1],
        ranges   = IRanges(start = bins[, 2], end = bins[, 3],
                           names = paste0(bins[, 1], ":", bins[, 2], "-", bins[, 3])),
        strand   = "+"
    )
    bamlist <- BamFileList(bam, yieldSize = 5e4)
    names(bamlist) <- basename(bam)
    counts <- summarizeOverlaps(features = coords, reads = bamlist,
                                ignore.strand = TRUE, singleEnd = single_end)
    dat <- as.data.frame(SummarizedExperiment::assay(counts))
    res <- data.frame(bin = rownames(dat), dat, check.names = FALSE)

    options(scipen = 99)
    total <- sum(as.numeric(dat[[1]]))
    if (total > min_reads) {
        write.table(res, out, sep = "\\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
        message(basename(bam), ": ", format(total, big.mark = ","), " reads counted in ",
                nrow(res), " bins")
    } else {
        # Bin column only: the merge step drops this sample, as CountRawReads would
        write.table(res[, "bin", drop = FALSE], out, sep = "\\t", quote = FALSE,
                    row.names = FALSE, col.names = TRUE)
        message("WARNING: ", basename(bam), " dropped: library size ",
                format(total, big.mark = ","), " <= 1e6 counted reads (as in ChIPseqSpikeInFree)")
    }
    END_R

    Rscript count_bam.R "${meta.id}.bam" "${chrom}" "${meta.id}" ${single_end} ${bin_size}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(paste(R.version\$major, R.version\$minor, sep = "."))' 2> /dev/null | tail -n 1)
        chipseqspikeinfree: \$(Rscript -e 'cat(as.character(packageVersion("ChIPseqSpikeInFree")))' 2> /dev/null | tail -n 1)
    END_VERSIONS
    """

    stub:
    """
    printf 'bin\\t%s.bam\\n' "${meta.id}" > ${meta.id}_rawCounts.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: stub
        chipseqspikeinfree: stub
    END_VERSIONS
    """
}
