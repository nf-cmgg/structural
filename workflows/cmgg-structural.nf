/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { fromSamplesheet; paramsSummaryMap } from 'plugin/nf-validation'

def summary_params = paramsSummaryMap(workflow)

// Check callers (see lib/GlobalVariables.groovy for the list of supported callers)
def callers = params.callers.tokenize(",")

for (caller in callers) {
    if(!(caller in GlobalVariables.allCallers)) { error("The caller '${caller}' is not supported please specify a comma delimited list with on or more of the following callers: ${GlobalVariables.allCallers}".toString()) }
}

def sv_callers_to_use = callers.intersect(GlobalVariables.svCallers)

if (sv_callers_to_use && params.callers_support > sv_callers_to_use.size()) {
    error("The --callers_support parameter (${params.callers_support}) is higher than the amount of SV callers in --callers (${sv_callers_to_use.size()}). Please adjust --callers_support to a value lower of equal to the amount of SV callers to use.")
}

if ("qdnaseq" in callers && !params.qdnaseq_reference) {
    error("Please give the qDNAseq reference using --qdnaseq_reference")
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { BAM_PREPARE_SAMTOOLS                  } from '../subworkflows/local/bam_prepare_samtools/main'
include { BAM_SV_CALLING                        } from '../subworkflows/local/bam_sv_calling/main'
include { BAM_CNV_CALLING                       } from '../subworkflows/local/bam_cnv_calling/main'
include { VCF_ANNOTATE_VEP_ANNOTSV_VCFANNO      } from '../subworkflows/local/vcf_annotate_vep_annotsv_vcfanno/main'
include { BAM_REPEAT_ESTIMATION_EXPANSIONHUNTER } from '../subworkflows/local/bam_repeat_estimation_expansionhunter/main'
include { VCF_CONCAT_BCFTOOLS                   } from '../subworkflows/local/vcf_concat_bcftools/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { TABIX_BGZIPTABIX                  } from '../modules/nf-core/tabix/bgziptabix/main'
include { SAMTOOLS_FAIDX                    } from '../modules/nf-core/samtools/faidx/main'
include { BWA_INDEX                         } from '../modules/nf-core/bwa/index/main'
include { ANNOTSV_INSTALLANNOTATIONS        } from '../modules/nf-core/annotsv/installannotations/main'
include { UNTAR as UNTAR_ANNOTSV            } from '../modules/nf-core/untar/main'
include { UNTAR as UNTAR_BWA                } from '../modules/nf-core/untar/main'
include { MULTIQC                           } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS       } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow CMGGSTRUCTURAL {

    ch_versions = Channel.empty()
    ch_reports  = Channel.empty()
    ch_outputs  = Channel.empty()
    count_types = 0

    //
    // Create input channels from parameters
    //

    ch_fasta_ready              = Channel.fromPath(params.fasta).map{[[id:'fasta'], it]}.collect()
    ch_fai                      = params.fai ?                      Channel.fromPath(params.fai).map{[[id:"fai"],it]}.collect() : null
    ch_bwa_index                = params.bwa ?                      Channel.fromPath(params.bwa).map{[[id:"bwa"],it]}.collect() : null
    ch_vep_cache                = params.vep_cache ?                Channel.fromPath(params.vep_cache).map{[[id:"vep_cache"],it]}.collect() : []
    ch_annotsv_annotations      = params.annotsv_annotations ?      Channel.fromPath(params.annotsv_annotations).map{[[id:"annotsv_annotations"], it]}.collect() :  null
    ch_annotsv_candidate_genes  = params.annotsv_candidate_genes ?  Channel.fromPath(params.annotsv_candidate_genes).map{[[], it]}.collect() : [[],[]]
    ch_annotsv_gene_transcripts = params.annotsv_gene_transcripts ? Channel.fromPath(params.annotsv_gene_transcripts).map{[[], it]}.collect() : [[],[]]
    ch_vcfanno_lua              = params.vcfanno_lua ?              Channel.fromPath(params.vcfanno_lua).collect() : []
    ch_catalog                  = params.expansionhunter_catalog ?  Channel.fromPath(params.expansionhunter_catalog).map{[[id:'catalog'], it]}.collect() : [[id:'catalog'],[file("https://github.com/Illumina/ExpansionHunter/raw/master/variant_catalog/grch38/variant_catalog.json", checkIfExists:true)]]    
    ch_qdnaseq_reference        = params.qdnaseq_reference ?        Channel.fromPath(params.qdnaseq_reference).map{[[id:'qdnaseq'], it]}.collect() : [[],[]]    
    ch_wisecondorx_reference    = params.wisecondorx_reference ?    Channel.fromPath(params.wisecondorx_reference).map{[[id:'wisecondorx'], it]}.collect() : [[],[]]    
    ch_blacklist                = params.blacklist ?                Channel.fromPath(params.blacklist).map{[[id:'blacklist'], it]}.collect() : [[],[]]    

    val_vcfanno_resources       = params.vcfanno_resources ?        params.vcfanno_resources.split(",").collect{file(it, checkIfExists:true)}.flatten() : []

    ch_vep_extra_files = []

    if(params.vep_phenotypes && params.phenotypes && params.phenotypes_tbi) {
        ch_vep_extra_files.add(file(params.phenotypes, checkIfExists:true))
        ch_vep_extra_files.add(file(params.phenotypes_tbi, checkIfExists:true))
    }
    else if(params.vep_phenotypes) {
        error("Please specify '--phenotypes PATH/TO/PHENOTYPES/FILE' and '--phenotypes_tbi PATH/TO/PHENOTYPES/INDEX/FILE' to use the Phenotypes VEP plugin.")
    }

    //
    // Create optional inputs
    //

    if(!ch_fai){
        SAMTOOLS_FAIDX(
            ch_fasta_ready
        )

        ch_versions  = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
        ch_fai_ready = SAMTOOLS_FAIDX.out.fai.map{[[id:'fai'], it]}.collect()
    }
    else {
        ch_fai_ready = ch_fai
    }

    if(!ch_bwa_index && "gridss" in callers){
        BWA_INDEX(
            ch_fasta_ready
        )

        ch_versions        = ch_versions.mix(BWA_INDEX.out.versions)
        ch_bwa_index_ready = BWA_INDEX.out.index.map{[[id:'bwa'], it[1]]}.collect()
    }
    else if(ch_bwa_index && params.bwa.endsWith(".tar.gz") && "gridss" in callers) {
        UNTAR_BWA(
            ch_bwa_index
        )
        ch_versions = ch_versions.mix(UNTAR_BWA.out.versions)

        UNTAR_BWA.out.untar
            .collect()
            .set { ch_bwa_index_ready }
    }
    else {
        ch_bwa_index_ready = ch_bwa_index
    }

    if(params.annotate && !ch_annotsv_annotations && callers.intersect(GlobalVariables.svCallers)) {
        ANNOTSV_INSTALLANNOTATIONS()
        ch_versions = ch_versions.mix(ANNOTSV_INSTALLANNOTATIONS.out.versions)

        ANNOTSV_INSTALLANNOTATIONS.out.annotations
            .map { [[id:"annotations"], it] }
            .collect()
            .set { ch_annotsv_annotations_ready }
    } 
    else if(params.annotate && callers.intersect(GlobalVariables.svCallers) && params.annotsv_annotations.endsWith(".tar.gz")) {
        UNTAR_ANNOTSV(
            ch_annotsv_annotations
        )
        ch_versions = ch_versions.mix(UNTAR_ANNOTSV.out.versions)

        UNTAR_ANNOTSV.out.untar
            .collect()
            .set { ch_annotsv_annotations_ready }
    }
    else {
        ch_annotsv_annotations_ready = ch_annotsv_annotations
    }

    //
    // Create the input channel
    //

    Channel.fromSamplesheet("input", immutable_meta:false)
        .multiMap({ meta, cram, crai, small_variants ->
            crams:          [ meta, cram, crai ]
            small_variants: [ meta, small_variants ]
        })
        .set { ch_inputs }

    //
    // Prepare the inputs
    //

    BAM_PREPARE_SAMTOOLS(
        ch_inputs.crams,
        ch_fasta_ready,
        ch_fai_ready
    )
    ch_versions = ch_versions.mix(BAM_PREPARE_SAMTOOLS.out.versions)

    //
    // Call the variants
    //

    if(callers.intersect(GlobalVariables.svCallers)){

        count_types++

        BAM_SV_CALLING(
            BAM_PREPARE_SAMTOOLS.out.crams,
            ch_fasta_ready,
            ch_fai_ready,
            ch_bwa_index_ready
        )

        ch_versions = ch_versions.mix(BAM_SV_CALLING.out.versions)
        ch_reports  = ch_reports.mix(BAM_SV_CALLING.out.reports)

        //
        // Annotate the variants
        //

        if(params.annotate) {
            VCF_ANNOTATE_VEP_ANNOTSV_VCFANNO(
                BAM_SV_CALLING.out.vcfs,
                ch_inputs.small_variants,
                ch_fasta_ready,
                ch_fai_ready,
                ch_annotsv_annotations_ready,
                ch_annotsv_candidate_genes,
                ch_annotsv_gene_transcripts,
                ch_vep_cache,
                ch_vep_extra_files,
                ch_vcfanno_lua,
                val_vcfanno_resources
            )

            ch_reports  = ch_reports.mix(VCF_ANNOTATE_VEP_ANNOTSV_VCFANNO.out.reports)
            ch_versions = ch_versions.mix(VCF_ANNOTATE_VEP_ANNOTSV_VCFANNO.out.versions)
            ch_outputs  = ch_outputs.mix(VCF_ANNOTATE_VEP_ANNOTSV_VCFANNO.out.vcfs)
        } else {
            ch_outputs  = ch_outputs.mix(BAM_SV_CALLING.out.vcfs)
        }
    }

    //
    // Copy number calling
    //

    if(callers.intersect(GlobalVariables.cnvCallers)){

        BAM_CNV_CALLING(
            BAM_PREPARE_SAMTOOLS.out.crams,
            ch_fasta_ready,
            ch_fai_ready,
            ch_qdnaseq_reference,
            ch_wisecondorx_reference,
            ch_blacklist
        )
        ch_versions = ch_versions.mix(BAM_CNV_CALLING.out.versions)
    }

    //
    // Estimate repeat sizes
    //

    if(callers.intersect(GlobalVariables.repeatsCallers)){

        count_types++

        BAM_REPEAT_ESTIMATION_EXPANSIONHUNTER(
            BAM_PREPARE_SAMTOOLS.out.crams,
            ch_fasta_ready,
            ch_fai_ready,
            ch_catalog
        )
        ch_versions = ch_versions.mix(BAM_REPEAT_ESTIMATION_EXPANSIONHUNTER.out.versions)
        ch_outputs  = ch_outputs.mix(BAM_REPEAT_ESTIMATION_EXPANSIONHUNTER.out.vcfs)

    }

    //
    // Concatenate the VCF files from different types of analysis
    //

    if(count_types > 1 && params.concat_output) {

        VCF_CONCAT_BCFTOOLS(
            ch_outputs,
            count_types
        )
        ch_versions = ch_versions.mix(VCF_CONCAT_BCFTOOLS.out.versions)

    }

    //
    // Dump the software versions
    //

    CUSTOM_DUMPSOFTWAREVERSIONS(
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    ch_versions_yaml = CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect()

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowNfCmggStructural.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowNfCmggStructural.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_reports)
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
