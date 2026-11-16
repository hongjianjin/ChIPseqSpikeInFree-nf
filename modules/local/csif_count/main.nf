process CSIF_COUNT {
    tag "${meta.id}"
    label 'process_medium'
    label 'process_long'

    container "${params.chipseqspikeinfree_container}"

    input:
    val  meta
    path sample_meta, stageAs: 'sample_meta.txt'
    path bam_map    , stageAs: 'bam_map.txt'
    path bams       , stageAs: 'bams/*'
    path chrom_sizes, stageAs: 'chrom_sizes/*'
    val  genome

    output:
    tuple val(meta), path("${prefix}_sample_meta.txt"), path("${prefix}_parsedMatrix.txt"), emit: parsed
    path "${prefix}_rawCounts.txt"                                         , emit: counts
    path "${prefix}_count.log"                                             , emit: log
    path "versions.yml"                                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    prefix         = task.ext.prefix ?: "${meta.id}"
    def genome_arg = chrom_sizes ? chrom_sizes : genome
    """
    cp sample_meta.txt ${prefix}_sample_meta.txt

    # Link each BAM as <sample>.bam: ChIPseqSpikeInFree uses the BAM file name as sample ID
    mkdir -p input
    : > bam.lst
    while IFS=\$'\\t' read -r sample bam; do
        [ -n "\$sample" ] || continue
        ln -s "../bams/\$bam" "input/\${sample}.bam"
        echo "\$PWD/input/\${sample}.bam" >> bam.lst
    done < bam_map.txt

    CSIF_SCRIPT=\$(Rscript -e 'cat(system.file("scripts/run_ChIPseqSpikeInFree.R", package = "ChIPseqSpikeInFree"))' 2> /dev/null | tail -n 1)
    if [ ! -f "\$CSIF_SCRIPT" ]; then
        echo "ERROR: run_ChIPseqSpikeInFree.R not found; is ChIPseqSpikeInFree installed in the container?" >&2
        exit 1
    fi

    # Standard mode counts reads, parses the counts, then calibrates with default cutoffs.
    # Only the counts and parsed matrix are kept here; CSIF_SCALE does the calibration, so a
    # calibration error at default cutoffs must not discard the (expensive) counts.
    status=0
    Rscript "\$CSIF_SCRIPT" \\
        --bamList bam.lst \\
        --genome ${genome_arg} \\
        --metaData ${prefix}_sample_meta.txt \\
        --prefix ${prefix} \\
        --ncores ${task.cpus} \\
        ${args} \\
        > ${prefix}_count.log 2>&1 || status=\$?
    cat ${prefix}_count.log

    if [ ! -s ${prefix}_rawCounts.txt ] || [ ! -s ${prefix}_parsedMatrix.txt ]; then
        echo "ERROR: read counting/parsing failed (exit \$status); see ${prefix}_count.log" >&2
        exit \$(( status == 0 ? 1 : status ))
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(paste(R.version\$major, R.version\$minor, sep = "."))' 2> /dev/null | tail -n 1)
        chipseqspikeinfree: \$(Rscript -e 'cat(as.character(packageVersion("ChIPseqSpikeInFree")))' 2> /dev/null | tail -n 1)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    cp sample_meta.txt ${prefix}_sample_meta.txt
    touch ${prefix}_rawCounts.txt ${prefix}_parsedMatrix.txt ${prefix}_count.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: stub
        chipseqspikeinfree: stub
    END_VERSIONS
    """
}
