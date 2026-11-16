/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CSIF_COUNT_BAM    } from '../modules/local/csif_count_bam/main'
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
    ch_samplesheet // channel: [ val(meta), path(bam) ]

    main:
    ch_versions    = Channel.empty()
    ch_chrom_sizes = params.chrom_sizes ? file(params.chrom_sizes, checkIfExists: true) : []

    //
    // Metadata in the layout ChIPseqSpikeInFree expects. ID is <sample>.bam, the name each
    // BAM is linked to before counting, so results are labelled by sample name. COLOR is
    // written only when every row has one.
    //
    ch_metadata = ch_samplesheet
        .map { meta, bam -> meta }
        .toSortedList { a, b -> a.id <=> b.id }
        .map { metas ->
            def with_color = metas.every { it.color }
            def header     = [ 'ID', 'ANTIBODY', 'GROUP' ] + (with_color ? [ 'COLOR' ] : [])
            def lines      = metas.collect { m ->
                ( [ "${m.id}.bam", m.antibody, m.group ] + (with_color ? [ m.color ] : []) ).join('\t')
            }
            ( [ header.join('\t') ] + lines ).join('\n') + '\n'
        }
        .collectFile(name: "${params.cohort}_sample_meta.txt")

    //
    // MODULE: count reads in genome-wide bins, one job per BAM (in parallel)
    //
    CSIF_COUNT_BAM (
        ch_samplesheet,
        ch_chrom_sizes,
        params.genome,
        params.bin_size,
        params.paired_end
    )
    ch_versions = ch_versions.mix(CSIF_COUNT_BAM.out.versions.first())

    //
    // MODULE: merge the per-BAM counts into one bin x sample table
    //
    CSIF_MERGE_COUNTS (
        [ id: params.cohort ],
        CSIF_COUNT_BAM.out.counts.map { meta, counts -> counts }.collect()
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
