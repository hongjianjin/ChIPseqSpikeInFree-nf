#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ChIPseqSpikeInFree-bigwig-nf
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Spike-in free ChIP-seq scaling factors with ChIPseqSpikeInFree, from bigWig coverage
    https://github.com/stjude/ChIPseqSpikeInFree
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CHIPSEQSPIKEINFREE      } from './workflows/chipseqspikeinfree'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_chipseqspikeinfree_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_chipseqspikeinfree_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CHIPSEQSPIKEINFREE_NF {

    take:
    samplesheet // channel: [ val(meta), val(bigwig: local path or URL) ]

    main:
    CHIPSEQSPIKEINFREE (
        samplesheet
    )

    emit:
    scaling_factors = CHIPSEQSPIKEINFREE.out.scaling_factors
    versions        = CHIPSEQSPIKEINFREE.out.versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.input
    )

    CHIPSEQSPIKEINFREE_NF (
        PIPELINE_INITIALISATION.out.samplesheet
    )

    PIPELINE_COMPLETION (
        params.outdir
    )
}
