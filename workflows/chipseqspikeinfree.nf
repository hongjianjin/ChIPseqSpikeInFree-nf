/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CSIF_BIN_BIGWIG   } from '../modules/local/csif_bin_bigwig/main'
include { CSIF_MERGE_COUNTS } from '../modules/local/csif_merge_counts/main'
include { CSIF_PARSE        } from '../modules/local/csif_parse/main'
include { CSIF_SCALE        } from '../modules/local/csif_scale/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CHIPSEQSPIKEINFREE {

    take:
    ch_samplesheet // channel: [ val(meta), val(bigwig: local path or URL string) ]

    main:
    ch_versions = Channel.empty()

    // Chromosome sizes always travel as a file: the bundled one for --genome, or --chrom_sizes
    ch_chrom_sizes = params.chrom_sizes
        ? file(params.chrom_sizes, checkIfExists: true)
        : file("${projectDir}/assets/chrom_sizes/${params.genome}.chrom.sizes", checkIfExists: true)

    //
    // URLs are passed as values so they are read in place with HTTP range requests; local
    // files are staged as paths.
    //
    ch_bigwig = ch_samplesheet.map { meta, bigwig ->
        def is_url = bigwig ==~ /^(https?|ftp):\/\/.+/
        [ meta, is_url ? bigwig : '', is_url ? [] : file(bigwig, checkIfExists: true) ]
    }

    //
    // Metadata in the layout ChIPseqSpikeInFree expects. ID is <sample>.bw, the column name
    // written by CSIF_BIN_BIGWIG, so results are labelled by sample name. COLOR is written only
    // when every row has one.
    //
    ch_metadata = ch_samplesheet
        .map { meta, bigwig -> meta }
        .toSortedList { a, b -> a.id <=> b.id }
        .map { metas ->
            def with_color = metas.every { it.color }
            def header     = [ 'ID', 'ANTIBODY', 'GROUP' ] + (with_color ? [ 'COLOR' ] : [])
            def lines      = metas.collect { m ->
                ( [ "${m.id}.bw", m.antibody, m.group ] + (with_color ? [ m.color ] : []) ).join('\t')
            }
            ( [ header.join('\t') ] + lines ).join('\n') + '\n'
        }
        .collectFile(name: "${params.cohort}_sample_meta.txt")

    //
    // MODULE: sum bigWig signal in genome-wide bins, one job per bigWig (in parallel)
    //
    CSIF_BIN_BIGWIG (
        ch_bigwig,
        ch_chrom_sizes,
        params.genome,
        params.bin_size,
        params.pseudo_depth,
        params.approx_summary
    )
    ch_versions = ch_versions.mix(CSIF_BIN_BIGWIG.out.versions.first())

    //
    // MODULE: merge the per-bigWig pseudo-counts into one bin x sample table
    //
    CSIF_MERGE_COUNTS (
        [ id: params.cohort ],
        CSIF_BIN_BIGWIG.out.counts.map { meta, counts -> counts }.collect()
    )
    ch_versions = ch_versions.mix(CSIF_MERGE_COUNTS.out.versions)

    //
    // MODULE: parse the merged counts into cumulative distributions
    //
    CSIF_PARSE (
        [ id: params.cohort ],
        ch_metadata,
        CSIF_MERGE_COUNTS.out.counts,
        ch_chrom_sizes,
        params.genome
    )
    ch_versions = ch_versions.mix(CSIF_PARSE.out.versions)

    //
    // MODULE: scaling factors, QC calls and plots from the parsed matrix
    // (re-runs alone when only calibration parameters change)
    //
    CSIF_SCALE (
        CSIF_PARSE.out.parsed
    )
    ch_versions = ch_versions.mix(CSIF_SCALE.out.versions)

    //
    // Collate software versions
    //
    ch_versions
        .unique()
        .collectFile(
            storeDir : "${params.outdir}/pipeline_info",
            name     : 'chipseqspikeinfree_software_versions.yml',
            sort     : true,
            newLine  : true
        )

    emit:
    scaling_factors = CSIF_SCALE.out.sf // channel: [ val(meta), path(SF.txt) ]
    versions        = ch_versions       // channel: [ path(versions.yml) ]
}
