#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-cmgg/structural
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-cmgg/structural
----------------------------------------------------------------------------------------
*/

include { validateParameters; paramsHelp; paramsSummaryLog } from 'plugin/nf-validation'

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.fasta                    = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.fai                      = WorkflowMain.getGenomeAttribute(params, 'fai')
params.vep_cache                = WorkflowMain.getGenomeAttribute(params, 'vep_cache')
params.bwa                      = WorkflowMain.getGenomeAttribute(params, 'bwa')
params.phenotypes               = WorkflowMain.getGenomeAttribute(params, 'phenotypes')
params.phenotypes_tbi           = WorkflowMain.getGenomeAttribute(params, 'phenotypes_tbi')
params.annotsv_annotations      = WorkflowMain.getGenomeAttribute(params, 'annotsv_annotations')
params.expansionhunter_catalog  = WorkflowMain.getGenomeAttribute(params, 'expansionhunter_catalog')
params.qdnaseq_male             = WorkflowMain.getGenomeAttribute(params, "qdnaseq_male_${params.qdnaseq_bin_size.toInteger() / 1000}kbp".toString())
params.qdnaseq_female           = WorkflowMain.getGenomeAttribute(params, "qdnaseq_female_${params.qdnaseq_bin_size.toInteger() / 1000}kbp".toString())
params.wisecondorx_reference    = WorkflowMain.getGenomeAttribute(params, 'wisecondorx_reference')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Print help message
if (params.help) {
    def String command = "nextflow run nf-cmgg/structural --input <input csv/tsv/yaml> --outdir <output folder>"
    log.info paramsHelp(command)
    exit 0
}

// Print parameter summary log to screen
log.info paramsSummaryLog(workflow)

// Validate input parameters
validateParameters()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CMGGSTRUCTURAL } from './workflows/cmgg-structural'

//
// WORKFLOW: Run main nf-cmgg/structural analysis pipeline
//
workflow CMGG_CMGGSTRUCTURAL {
    CMGGSTRUCTURAL ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    CMGG_CMGGSTRUCTURAL ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
