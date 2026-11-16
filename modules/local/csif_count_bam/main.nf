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
    suppressMessages(library(ChIPseqSpikeInFree))
    args <- commandArgs(trailingOnly = TRUE)
    bam <- args[1]; chrom <- args[2]; prefix <- args[3]
    single_end <- as.logical(args[4]); bin_size <- as.numeric(args[5])
    out <- paste0(prefix, "_rawCounts.txt")

    # Same counting as standard mode of run_ChIPseqSpikeInFree.R, one BAM at a time
    res <- tryCatch(
        CountRawReads(bamFiles = bam, chromFile = chrom, prefix = prefix,
                      singleEnd = single_end, binSize = bin_size),
        error = function(e) e
    )
    if (inherits(res, "error")) {
        # CountRawReads drops BAMs with < 1e6 counted reads; it then writes the bin column only
        # and fails when returning. Keep that file so the merge step drops the sample, as the
        # single-job run would; any other failure is fatal.
        if (file.exists(out)) {
            message("WARNING: ", basename(bam), " dropped by ChIPseqSpikeInFree ",
                    "(library size < 1e6 reads or no counts): ", conditionMessage(res))
        } else {
            stop(res)
        }
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
