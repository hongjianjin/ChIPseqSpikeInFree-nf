process CSIF_BIN_BIGWIG {
    tag "${meta.id}"
    label 'process_single'
    label 'error_retry'

    container "${params.chipseqspikeinfree_container}"

    input:
    // Exactly one of bigwig_url / bigwig_file is set: http(s) and ftp URLs are read in place with
    // HTTP range requests (not downloaded); anything else is a staged local file.
    tuple val(meta), val(bigwig_url), path(bigwig_file, stageAs: 'bigwig/*')
    path  chrom_sizes, stageAs: 'chrom_sizes/*'
    val   genome
    val   bin_size
    val   pseudo_depth
    val   approx

    output:
    tuple val(meta), path("${meta.id}_rawCounts.txt"), emit: counts
    path "${meta.id}_bigwig_qc.json"                 , emit: qc
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def source     = bigwig_url ? bigwig_url : bigwig_file
    def chrom      = chrom_sizes ? chrom_sizes : genome
    def reads      = meta.reads ? "--reads ${meta.reads}" : "--pseudo-depth ${pseudo_depth}"
    def approx_arg = approx ? '--approx' : ''
    """
    # Same genome-wide bins as ChIPseqSpikeInFree::GenerateBins(); the signal mass of each bin is
    # rescaled to pseudo read counts (ParseReadCounts() divides by the sample total, so the
    # bigWig's own scaling cancels). The column is named <sample>.bw, the ID used in the metadata.
    csif_bigwig_bin.py \\
        --bigwig "${source}" \\
        --id "${meta.id}.bw" \\
        --genome "${chrom}" \\
        --prefix "${meta.id}" \\
        --bin-size ${bin_size} \\
        ${reads} \\
        ${approx_arg} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
        pybigtools: \$(python3 -c 'import importlib.metadata as m; print(m.version("pybigtools"))')
    END_VERSIONS
    """

    stub:
    """
    printf 'bin\\t%s.bw\\n' "${meta.id}" > ${meta.id}_rawCounts.txt
    echo '{}' > ${meta.id}_bigwig_qc.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
        pybigtools: stub
    END_VERSIONS
    """
}
