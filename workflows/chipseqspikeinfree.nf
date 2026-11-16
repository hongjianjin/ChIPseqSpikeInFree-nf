/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CSIF_COUNT } from '../modules/local/csif_count/main'
include { CSIF_SCALE } from '../modules/local/csif_scale/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CHIPSEQSPIKEINFREE {

    take:
    ch_samplesheet // channel: [ val(meta), path(bam) ]

    main:
    ch_versions = Channel.empty()

    //
    // All BAMs of the cohort go into one run. Metadata in the layout ChIPseqSpikeInFree
    // expects: ID is the file name CSIF_COUNT links each BAM to (<sample>.bam), so
    // results are labelled by sample name. COLOR is written only when every row has one.
    //
    ch_rows = ch_samplesheet
        .map { meta, bam -> [ meta, bam ] }
        .toList()
        .map { rows ->
            def dups = rows.groupBy { it[1].name }.findAll { k, v -> v.size() > 1 }
            if (dups) {
                error("BAM file names must be unique across samples; repeated: " +
                    dups.collect { k, v -> "${k} (${v.collect { it[0].id }.join(', ')})" }.join('; '))
            }
            rows.sort { it[0].id }
        }

    ch_metadata = ch_rows
        .map { rows ->
            def with_color = rows.every { it[0].color }
            def header     = [ 'ID', 'ANTIBODY', 'GROUP' ] + (with_color ? [ 'COLOR' ] : [])
            def lines      = rows.collect { meta, bam ->
                ( [ "${meta.id}.bam", meta.antibody, meta.group ] + (with_color ? [ meta.color ] : []) ).join('\t')
            }
            ( [ header.join('\t') ] + lines ).join('\n') + '\n'
        }
        .collectFile(name: "${params.cohort}_sample_meta.txt")

    // sample -> staged BAM file name, used to link each BAM as <sample>.bam
    ch_bam_map = ch_rows
        .map { rows -> rows.collect { meta, bam -> "${meta.id}\t${bam.name}" }.join('\n') + '\n' }
        .collectFile(name: "${params.cohort}_bam_map.txt")

    ch_bams = ch_rows.map { rows -> rows.collect { it[1] } }

    ch_chrom_sizes = params.chrom_sizes ? file(params.chrom_sizes, checkIfExists: true) : []

    //
    // MODULE: count reads in genome-wide bins and parse the counts
    //
    CSIF_COUNT (
        [ id: params.cohort ],
        ch_metadata,
        ch_bam_map,
        ch_bams,
        ch_chrom_sizes,
        params.genome
    )
    ch_versions = ch_versions.mix(CSIF_COUNT.out.versions)

    //
    // MODULE: scaling factors, QC calls and plots from the parsed matrix
    // (re-runs alone when only calibration parameters change)
    //
    CSIF_SCALE (
        CSIF_COUNT.out.parsed
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
