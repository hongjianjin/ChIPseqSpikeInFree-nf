process CSIF_PARSE {
    tag "${meta.id}"
    label 'process_medium'

    container "${params.chipseqspikeinfree_container}"

    input:
    val  meta
    path sample_meta, stageAs: 'sample_meta.txt'
    path raw_counts , stageAs: 'merged/*'
    path chrom_sizes, stageAs: 'chrom_sizes/*'
    val  genome

    output:
    tuple val(meta), path("${prefix}_sample_meta.txt"), path("${prefix}_parsedMatrix.txt"), emit: parsed
    path "${prefix}_parse.log"                                                             , emit: log
    path "versions.yml"                                                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    prefix         = task.ext.prefix ?: "${meta.id}"
    def genome_arg = chrom_sizes ? chrom_sizes : genome
    """
    cp sample_meta.txt ${prefix}_sample_meta.txt

    # run_ChIPseqSpikeInFree.R loads an existing <prefix>_rawCounts.txt instead of counting,
    # so this step only parses the merged per-BAM counts. It still requires a BAM list file.
    ln -s "${raw_counts}" ${prefix}_rawCounts.txt
    tail -n +2 ${prefix}_sample_meta.txt | cut -f 1 > bam.lst

    CSIF_SCRIPT=\$(Rscript -e 'cat(system.file("scripts/run_ChIPseqSpikeInFree.R", package = "ChIPseqSpikeInFree"))' 2> /dev/null | tail -n 1)
    if [ ! -f "\$CSIF_SCRIPT" ]; then
        echo "ERROR: run_ChIPseqSpikeInFree.R not found; is ChIPseqSpikeInFree installed in the container?" >&2
        exit 1
    fi

    # Standard mode parses the counts, then calibrates with default cutoffs. Only the parsed
    # matrix is kept; CSIF_SCALE does the calibration, so a calibration error at default
    # cutoffs must not discard the parsed matrix.
    status=0
    Rscript "\$CSIF_SCRIPT" \\
        --bamList bam.lst \\
        --genome ${genome_arg} \\
        --metaData ${prefix}_sample_meta.txt \\
        --prefix ${prefix} \\
        --ncores ${task.cpus} \\
        ${args} \\
        > ${prefix}_parse.log 2>&1 || status=\$?
    cat ${prefix}_parse.log

    if [ ! -s ${prefix}_parsedMatrix.txt ]; then
        echo "ERROR: parsing raw counts failed (exit \$status); see ${prefix}_parse.log" >&2
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
    touch ${prefix}_parsedMatrix.txt ${prefix}_parse.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: stub
        chipseqspikeinfree: stub
    END_VERSIONS
    """
}
