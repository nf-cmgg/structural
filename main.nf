#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-cmgg/structural
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-cmgg/structural
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { STRUCTURAL              } from './workflows/structural'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_structural_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_structural_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_structural_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.fasta                    = getGenomeAttribute('fasta')
params.fai                      = getGenomeAttribute('fai')
params.vep_cache                = getGenomeAttribute('vep_cache')
// params.bwa                      = getGenomeAttribute('bwa')
params.annotsv_annotations      = getGenomeAttribute('annotsv_annotations')
params.expansionhunter_catalog  = getGenomeAttribute('expansionhunter_catalog')
params.qdnaseq_male             = getGenomeAttribute("qdnaseq_male_${params.qdnaseq_bin_size.toInteger() / 1000}kbp".toString())
params.qdnaseq_female           = getGenomeAttribute("qdnaseq_female_${params.qdnaseq_bin_size.toInteger() / 1000}kbp".toString())
params.wisecondorx_reference    = getGenomeAttribute('wisecondorx_reference')


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    def ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    def ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
    def ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
    def ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

    //
    // WORKFLOW: Run main workflow
    //
    STRUCTURAL (
        // channels
        PIPELINE_INITIALISATION.out.samplesheet,
        ch_multiqc_config,
        ch_multiqc_custom_config,
        ch_multiqc_logo,
        ch_multiqc_custom_methods_description,

        // files
        params.fasta,
        params.fai,
        params.expansionhunter_catalog ?: "https://github.com/Illumina/ExpansionHunter/raw/master/variant_catalog/grch38/variant_catalog.json",
        params.qdnaseq_female,
        params.qdnaseq_male,
        params.wisecondorx_reference,
        params.vep_cache,
        params.annotsv_annotations,
        params.annotsv_candidate_genes,
        params.annotsv_gene_transcripts,
        params.vcfanno_lua,
        params.vcfanno_resources,
        params.vcfanno_toml,
        params.blacklist,
        params.manta_config ?: "${projectDir}/assets/manta_config.ini",
        "${projectDir}/assets/svync",
        "${projectDir}/assets/bedgovcf",
        "${projectDir}/assets/vcfanno",

        // booleans
        params.annotate,
        params.concat_output,

        // values
        params.callers,
        params.sv_callers_support,
        params.cnv_callers_support,
        params.genome,
        params.species,
        params.vep_assembly,
        params.vep_cache_version,
        params.annotations_filter,
        params.outdir
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
        STRUCTURAL.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
