//
// Subworkflows with functionality specific to the ChIPseqSpikeInFree-nf pipeline
//

include { validateParameters } from 'plugin/nf-schema'
include { paramsSummaryLog   } from 'plugin/nf-schema'
include { samplesheetToList  } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version         // boolean: Display version and exit
    validate_params // boolean: Validate parameters against the schema at runtime
    input           // string:  Path to input samplesheet

    main:

    //
    // Print version and exit if required
    //
    if (version) {
        log.info "${workflow.manifest.name} ${getWorkflowVersion()}"
        System.exit(0)
    }

    //
    // Validate parameters and print summary to screen
    //
    if (validate_params) {
        validateParameters()
    }
    log.info paramsSummaryLog(workflow)

    //
    // Create channel from input file provided through params.input
    //
    Channel
        .fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet // channel: [ val(meta), path(bam) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    outdir // path: Path to output directory where results will be published

    main:
    // Keep `workflow` references out of the handler closure: Nextflow runs it with the
    // workflow metadata as delegate, so `workflow.success` there does not resolve.
    workflow.onComplete {
        completionSummary(outdir)
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to the troubleshooting docs in docs/usage.md"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Print a one-line completion summary
//
def completionSummary(outdir) {
    if (workflow.success) {
        log.info "-[${workflow.manifest.name}] Pipeline completed successfully. Results: ${outdir}"
    } else {
        log.error "-[${workflow.manifest.name}] Pipeline completed with errors"
    }
}

//
// Version string for the pipeline
//
def getWorkflowVersion() {
    def version_string = "" as String
    if (workflow.manifest.version) {
        def prefix_v = workflow.manifest.version[0] != 'v' ? 'v' : ''
        version_string += "${prefix_v}${workflow.manifest.version}"
    }
    if (workflow.commitId) {
        def git_shortsha = workflow.commitId.substring(0, 7)
        version_string += "-g${git_shortsha}"
    }
    return version_string
}
