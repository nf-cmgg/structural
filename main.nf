#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-cmgg/structural
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-cmgg/structural
----------------------------------------------------------------------------------------
*/

include { validateParameters; paramsHelp; paramsSummaryLog } from 'plugin/nf-schema'

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_structural_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { STRUCTURAL              } from './workflows/structural'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_structural_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_structural_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCMGG_STRUCTURAL {

    take:
    samplesheet     // channel: samplesheet read in from --input
    pipeline_params //     map: The parameters needed for the pipeline to run

    main:

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        CONFIG FILES
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config   = pipeline_params.multiqc_config ? Channel.fromPath( pipeline_params.multiqc_config, checkIfExists: true ) : Channel.empty()
    ch_multiqc_logo            = pipeline_params.multiqc_logo   ? Channel.fromPath( pipeline_params.multiqc_logo, checkIfExists: true ) : Channel.empty()
    ch_multiqc_custom_methods_description = pipeline_params.multiqc_methods_description ? file(pipeline_params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)


    //
    // WORKFLOW: Run pipeline
    //
    STRUCTURAL (
        // channels
        samplesheet,
        ch_multiqc_config,
        ch_multiqc_custom_config,
        ch_multiqc_logo,
        ch_multiqc_custom_methods_description,

        // file inputs
        pipeline_params.fasta,
        pipeline_params.fai,
        pipeline_params.expansionhunter_catalog ?: "https://github.com/Illumina/ExpansionHunter/raw/master/variant_catalog/grch38/variant_catalog.json",
        pipeline_params.qdnaseq_female,
        pipeline_params.qdnaseq_male,
        pipeline_params.wisecondorx_reference,
        pipeline_params.vep_cache,
        pipeline_params.annotsv_annotations,
        pipeline_params.annotsv_candidate_genes,
        pipeline_params.annotsv_gene_transcripts,
        pipeline_params.vcfanno_lua,
        pipeline_params.vcfanno_resources,
        pipeline_params.vcfanno_toml,
        pipeline_params.blacklist,
        pipeline_params.manta_config ?: "${projectDir}/assets/manta_config.ini",
        "${projectDir}/assets/svync",
        "${projectDir}/assets/bedgovcf",
        "${projectDir}/assets/vcfanno",

        // boolean inputs
        pipeline_params.annotate,
        pipeline_params.concat_output,

        // value inputs
        pipeline_params.callers,
        pipeline_params.sv_callers_support,
        pipeline_params.cnv_callers_support,
        pipeline_params.genome,
        pipeline_params.species,
        pipeline_params.vep_assembly,
        pipeline_params.vep_cache_version,
        pipeline_params.annotations_filter,
        pipeline_params.outdir
    )

    emit:
    multiqc_report = STRUCTURAL.out.multiqc_report // channel: /path/to/multiqc_report.html

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    //
    // Initialize genome parameters
    //
    params.fasta                    = getGenomeAttribute('fasta', params.genome, params.genomes)
    params.fai                      = getGenomeAttribute('fai', params.genome, params.genomes)
    params.vep_cache                = getGenomeAttribute('vep_cache', params.genome, params.genomes)
    // params.bwa                      = getGenomeAttribute('bwa', params.genome, params.genomes)
    params.annotsv_annotations      = getGenomeAttribute('annotsv_annotations', params.genome, params.genomes)
    params.expansionhunter_catalog  = getGenomeAttribute('expansionhunter_catalog', params.genome, params.genomes)
    params.qdnaseq_male             = getGenomeAttribute("qdnaseq_male_${params.qdnaseq_bin_size.toInteger() / 1000}kbp".toString(), params.genome, params.genomes)
    params.qdnaseq_female           = getGenomeAttribute("qdnaseq_female_${params.qdnaseq_bin_size.toInteger() / 1000}kbp".toString(), params.genome, params.genomes)
    params.wisecondorx_reference    = getGenomeAttribute('wisecondorx_reference', params.genome, params.genomes)

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.genome,
        params.genomes
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCMGG_STRUCTURAL (
        PIPELINE_INITIALISATION.out.samplesheet,
        params
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCMGG_STRUCTURAL.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
