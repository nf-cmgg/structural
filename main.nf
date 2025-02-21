#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-cmgg/structural
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-cmgg/structural
----------------------------------------------------------------------------------------
*/

nextflow.preview.output = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { STRUCTURAL  } from './workflows/structural'
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
        params.bedpe,

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

    publish:
    STRUCTURAL.out.caller_vcfs      >> 'caller_vcfs'
    STRUCTURAL.out.sample_vcfs      >> 'sample_vcfs'
    STRUCTURAL.out.family_vcfs      >> 'family_vcfs'
    STRUCTURAL.out.qdnaseq_out      >> 'qdnaseq_out'
    STRUCTURAL.out.wisecondorx_out  >> 'wisecondorx_out'
    STRUCTURAL.out.bedpe            >> 'bedpe'
    STRUCTURAL.out.multiqc_report   >> 'multiqc'
    STRUCTURAL.out.multiqc_data     >> 'multiqc_data'

}

output {
    'caller_vcfs' {
        enabled params.output_callers
        path { meta, vcf, _tbi -> { file ->
            if(file == vcf.name) {
                return "${meta.sample}/${meta.caller}/${meta.sample}.vcf.gz"
            }
            return "${meta.sample}/${meta.caller}/${meta.sample}.vcf.gz.tbi"
        } }
    }
    'sample_vcfs' {
        path { meta, vcf, _tbi -> { file ->
            def base = "${meta.id}/${meta.id}${meta.variant_type ? '.' + meta.variant_type : ''}"
            if(file == vcf.name) {
                return "${base}.vcf.gz"
            }
            return "${base}.vcf.gz.tbi"
        } }
    }
    'family_vcfs' {
        path { meta, vcf, _tbi -> { file ->
            def base = "${meta.id}/${meta.id}${meta.variant_type ? '.' + meta.variant_type : ''}"
            if(file == vcf.name) {
                return "${base}.vcf.gz"
            }
            return "${base}.vcf.gz.tbi"
        } }
    }
    'qdnaseq_out' {
        path { meta, _bed -> { file ->
            if(file == "statistics.out") {
                return "${meta.id}/${meta.id}.qdnaseq.statistics.out"
            }
            def new_name = file.replaceFirst(meta.id, "${meta.id}.qdnaseq")
            return "${meta.id}/${new_name}"
        } }
    }
    'wisecondorx_out' {
        path { meta, _bed -> { file ->
            if(file.endsWith(".png")) {
                return "${meta.id}/${meta.id}.wisecondorx.${file}"
            }
            def new_name = file.replaceFirst(meta.id, "${meta.id}.wisecondorx")
            return "${meta.id}/${new_name}"
        } }
    }
    'bedpe' {
        path { meta, _bedpe -> { file ->
            def base = "${meta.id}/${meta.id}${meta.variant_type ? '.' + meta.variant_type : ''}"
            return "${base}.bedpe"
        } }
    }
    'multiqc' {
        path { _report -> { _file -> "multiqc/multiqc_report.html"}}
    }
    'multiqc_data' {
        path { _folder -> { _file -> "multiqc/multiqc_data"}}
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
