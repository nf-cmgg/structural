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
params.dict                     = getGenomeAttribute('dict')
params.gtf                      = getGenomeAttribute('gtf')
params.vep_cache                = getGenomeAttribute('vep_cache')
// params.bwa                      = getGenomeAttribute('bwa')
// params.annotsv_annotations      = getGenomeAttribute('annotsv_annotations')
params.expansionhunter_catalog  = getGenomeAttribute('expansionhunter_catalog')
params.qdnaseq_male             = getGenomeAttribute("qdnaseq_male_${params.qdnaseq_bin_size.toInteger() / 1000}kbp".toString())
params.qdnaseq_female           = getGenomeAttribute("qdnaseq_female_${params.qdnaseq_bin_size.toInteger() / 1000}kbp".toString())
params.wisecondorx_reference    = getGenomeAttribute('wisecondorx_reference')
params.strvctvre_phylop         = getGenomeAttribute('strvctvre_phylop')
params.strvctvre_data           = getGenomeAttribute('strvctvre_data')


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
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    def ch_multiqc_config          = channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    def ch_multiqc_custom_config   = params.multiqc_config ? channel.fromPath( params.multiqc_config, checkIfExists: true ) : channel.empty()
    def ch_multiqc_logo            = params.multiqc_logo   ? channel.fromPath( params.multiqc_logo, checkIfExists: true ) : channel.empty()
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
        params.dict,
        params.gtf,
        params.expansionhunter_catalog ?: "https://github.com/Illumina/ExpansionHunter/raw/master/variant_catalog/grch38/variant_catalog.json",
        params.qdnaseq_female,
        params.qdnaseq_male,
        params.wisecondorx_reference,
        params.vep_cache,
        // params.annotsv_annotations,
        // params.annotsv_candidate_genes,
        // params.annotsv_gene_transcripts,
        params.vcfanno_lua,
        params.vcfanno_resources,
        params.vcfanno_toml,
        params.blacklist,
        params.manta_config ?: "${projectDir}/assets/manta_config.ini",
        "${projectDir}/assets/svync",
        "${projectDir}/assets/bedgovcf",
        "${projectDir}/assets/vcfanno",
        params.strvctvre_phylop,
        params.strvctvre_data,

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
        params.filter,
        params.outdir,
        params.annotate_tools ? params.annotate_tools.tokenize(",") : []
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
    caller_vcfs     = STRUCTURAL.out.caller_vcfs
    sample_vcfs     = STRUCTURAL.out.sample_vcfs
    family_vcfs     = STRUCTURAL.out.family_vcfs
    qdnaseq_out     = STRUCTURAL.out.qdnaseq_out
    wisecondorx_out = STRUCTURAL.out.wisecondorx_out
    bedpe           = STRUCTURAL.out.bedpe
    multiqc         = STRUCTURAL.out.multiqc_report
    multiqc_data    = STRUCTURAL.out.multiqc_data

}

output {
    caller_vcfs {
        enabled params.output_callers
        path { meta, vcf, tbi ->
            vcf >> "${meta.sample}/${meta.caller}/${meta.sample}.vcf.gz"
            tbi >> "${meta.sample}/${meta.caller}/${meta.sample}.vcf.gz.tbi"
        }
    }
    sample_vcfs {
        path { meta, vcf, tbi ->
            def base = "${meta.id}/${meta.id}${meta.variant_type ? '.' + meta.variant_type : ''}"
            vcf >> "${base}.vcf.gz"
            tbi >> "${base}.vcf.gz.tbi"
        }
    }
    family_vcfs {
        path { meta, vcf, tbi ->
            def base = "${meta.id}/${meta.id}${meta.variant_type ? '.' + meta.variant_type : ''}"
            vcf >> "${base}.vcf.gz"
            tbi >> "${base}.vcf.gz.tbi"
        }
    }
    qdnaseq_out {
        path { meta, _bed_qdnaseq -> "$meta.id/qdnaseq/"
            // def base_suffix = bed_qdnaseq.name.replace(meta.id, "${meta.id}.qdnaseq")
            // bed_qdnaseq >> bed_qdnaseq.name == "statistics.out" ?
            //     "${meta.id}/${meta.id}.qdnaseq.statistics.out" :
            //     "${meta.id}/${base_suffix}"
        }
    }
    wisecondorx_out {
        path { meta, _bed -> "$meta.id/wisecondorx/"
            // if(bed.name.endsWith(".png")) {
            //     bed >> "${meta.id}/${meta.id}.wisecondorx.${bed.name}"
            // } else {
            //     def new_name = bed.name.replaceFirst(meta.id, "${meta.id}.wisecondorx")
            //     bed >> "${meta.id}/${new_name}"
            // }
        }
    }
    bedpe {
        path { meta, bedpe ->
            def base = "${meta.id}/${meta.id}${meta.variant_type ? '.' + meta.variant_type : ''}"
            bedpe >> "${base}.bedpe"
        }
    }
    multiqc {
        path { "multiqc/" }
    }
    multiqc_data {
        path { "multiqc/" }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
