process CSIF_SCALE {
    tag "${meta.id}"
    label 'process_single'

    container "${params.chipseqspikeinfree_container}"

    input:
    tuple val(meta), path(sample_meta), path(parsed)

    output:
    tuple val(meta), path("${prefix}_SF.txt"), emit: sf
    path "${prefix}_distribution.pdf"        , emit: distribution, optional: true
    path "${prefix}_boxplot.pdf"             , emit: boxplot     , optional: true
    path "${prefix}_scale.log"               , emit: log
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    CSIF_SCRIPT=\$(Rscript -e 'cat(system.file("scripts/run_ChIPseqSpikeInFree.R", package = "ChIPseqSpikeInFree"))' 2> /dev/null | tail -n 1)
    if [ ! -f "\$CSIF_SCRIPT" ]; then
        echo "ERROR: run_ChIPseqSpikeInFree.R not found; is ChIPseqSpikeInFree installed in the container?" >&2
        exit 1
    fi

    # The parsed matrix must be named <prefix>_parsedMatrix.txt for replotting mode
    if [ "${parsed}" != "${meta.id}_parsedMatrix.txt" ]; then
        ln -sf "${parsed}" "${meta.id}_parsedMatrix.txt"
    fi

    # Replotting mode: scaling factors, QC calls and plots from the parsed matrix
    Rscript "\$CSIF_SCRIPT" \\
        --metaData ${sample_meta} \\
        --prefix ${meta.id} \\
        --replot TRUE \\
        ${args} \\
        2>&1 | tee ${prefix}_scale.log

    # Drop the "_replot" tag the script adds to its output names
    for f in ${meta.id}_replot_*; do
        [ -e "\$f" ] || continue
        mv "\$f" "${prefix}_\${f#${meta.id}_replot_}"
    done

    if [ ! -s ${prefix}_SF.txt ]; then
        echo "ERROR: ${prefix}_SF.txt was not produced; see ${prefix}_scale.log" >&2
        exit 1
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
    touch ${prefix}_SF.txt ${prefix}_distribution.pdf ${prefix}_boxplot.pdf ${prefix}_scale.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: stub
        chipseqspikeinfree: stub
    END_VERSIONS
    """
}
