/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowNfCmggStructural.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input,
    params.multiqc_config,
    params.fasta,
    params.fai,
    params.vep_cache,
    params.phenotypes,
    params.phenotypes_tbi,
    params.annotsv_annotations,
    params.vcfanno_toml,
    params.vcfanno_lua
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input, checkIfExists: true) } else { exit 1, 'Input samplesheet not specified!' }

// Check callers
def availableCallers = [
    "delly",
    // "whamg",
    "manta",
    // "gridss",
    "smoove"
]

for (caller in params.callers.tokenize(",")) {
    if(!(caller in availableCallers)) { error("The caller '${caller}' is not supported please specify a comma delimited list with on or more of the following callers: ${availableCallers}".toString()) }
}

if ("whamg" in params.callers.tokenize(",")) {
    error("Whamg currently isn't functional. This will be fixed in a further build of the pipeline")
}

if ("gridss" in params.callers.tokenize(",")) {
    error("Gridss currently isn't functional. This will be fixed in a further build of the pipeline")
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
include { BAM_STRUCTURAL_VARIANT_CALLING    } from '../subworkflows/local/bam_structural_variant_calling/main'
include { VCF_GENOTYPE_SV_PARAGRAPH         } from '../subworkflows/local/vcf_genotype_sv_paragraph/main'
include { VCF_ANNOTATE_VEP_ANNOTSV_VCFANNO  } from '../subworkflows/local/vcf_annotate_vep_annotsv_vcfanno/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { TABIX_BGZIPTABIX                  } from '../modules/nf-core/tabix/bgziptabix/main'
include { BEDTOOLS_SORT                     } from '../modules/nf-core/bedtools/sort/main'
include { SAMTOOLS_FAIDX                    } from '../modules/nf-core/samtools/faidx/main'
include { BWA_INDEX                         } from '../modules/nf-core/bwa/index/main'
include { ANNOTSV_INSTALLANNOTATIONS        } from '../modules/nf-core/annotsv/installannotations/main'
include { UNTAR                             } from '../modules/nf-core/untar/main'
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

    //
    // Create input channels from parameters
    //

    ch_fasta                    = Channel.fromPath(params.fasta).collect()
    ch_fai                      = params.fai ?                      Channel.fromPath(params.fai).collect() :                                                null
    ch_bwa_index                = params.bwa ?                      Channel.fromPath(params.bwa).map {[[],it]}.collect() :                                  null
    ch_vep_cache                = params.vep_cache ?                Channel.fromPath(params.vep_cache).collect() :                                          []
    ch_annotsv_annotations      = params.annotsv_annotations ?      Channel.fromPath(params.annotsv_annotations).map{[[id:"annotations"], it]}.collect() :  null
    ch_annotsv_candidate_genes  = params.annotsv_candidate_genes ?  Channel.fromPath(params.annotsv_candidate_genes).map{[[], it]}.collect() :              [[],[]]
    ch_annotsv_gene_transcripts = params.annotsv_gene_transcripts ? Channel.fromPath(params.annotsv_gene_transcripts).map{[[], it]}.collect() :             [[],[]]
    ch_vcfanno_lua              = params.vcfanno_lua ?              Channel.fromPath(params.vcfanno_lua).collect() :                                        []
    val_vcfanno_resources       = params.vcfanno_resources ?        params.vcfanno_resources.split(",").collect{file(it, checkIfExists:true)}.flatten() :   []

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
            ch_fasta.map {[[],it]}
        )

        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
        ch_fai   = SAMTOOLS_FAIDX.out.fai.map { it[1] }.collect()
    }

    if(!ch_bwa_index && params.callers.contains("gridss")){
        BWA_INDEX(
            ch_fasta.map {[[id:'bwa'],it]}
        )

        ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
        ch_bwa_index = BWA_INDEX.out.index.collect()
    }

    if(params.annotate && !ch_annotsv_annotations) {
        ANNOTSV_INSTALLANNOTATIONS()
        ch_versions = ch_versions.mix(ANNOTSV_INSTALLANNOTATIONS.out.versions)

        ANNOTSV_INSTALLANNOTATIONS.out.annotations
            .map { [[id:"annotations"], it] }
            .collect()
            .set { ch_annotsv_annotations_ready }
    } 
    else if(params.annotate && params.annotsv_annotations.endsWith(".tar.gz")) {
        UNTAR(
            ch_annotsv_annotations
        )
        ch_versions = ch_versions.mix(UNTAR.out.versions)

        UNTAR.out.untar
            .collect()
            .set { ch_annotsv_annotations_ready }
    }
    else {
        ch_annotsv_annotations_ready = ch_annotsv_annotations
    }

    //
    // Create the input channel
    //

    SamplesheetConversion.convert(ch_input, file("${projectDir}/assets/schema_input.json"))
        .map { meta, cram, crai, bed, ped, small_variants ->
            new_meta = meta + [family:meta.family ?: meta.id]
            [ new_meta, cram, crai, bed, ped, small_variants ]
        }
        .tap { original_samplesheet }
        .map { meta, cram, crai, bed, ped, small_variants ->
            [ meta.family, 1 ]
        }
        .groupTuple()
        .map { family, one ->
            [ family, one.sum() ]
        }
        .combine(
            original_samplesheet.map {
                [ it[0].family ] + it
            },
            by:0
        )
        .multiMap({ family, family_count, meta, cram, crai, bed, ped, small_variants ->
            new_meta = meta + [family_count:family_count]
            family_meta = [
                id: meta.family,
                family: meta.family,
                family_count: family_count
            ]
            bed: [ new_meta, bed ]
            crams: [ new_meta, cram, crai ]
            small_variants: [ family_meta, small_variants ]
        })
        .set { ch_inputs }

    //
    // Use one small variants file per family
    //

    ch_inputs.small_variants
        .groupTuple() // No size needed here because no process has been run with small variant VCF files before this
        .map { meta, vcfs ->
            // Find the first VCF file and return that one for the family ([] if no VCF is given for the family)
            [ meta, vcfs.find { it != [] } ?: [] ]
        }
        .set { ch_small_variants_ready }

    //
    // Prepare the BED files
    //

    ch_inputs.bed
        .branch { meta, bed ->
            bed: bed
            no_bed: !bed
                return [ meta, [], [], [] ]
        }
        .set { ch_all_beds }

    BEDTOOLS_SORT(
        ch_all_beds.bed,
        []
    )

    ch_versions = ch_versions.mix(BEDTOOLS_SORT.out.versions)

    TABIX_BGZIPTABIX(
        BEDTOOLS_SORT.out.sorted
    )
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions)

    BEDTOOLS_SORT.out.sorted
        .join(TABIX_BGZIPTABIX.out.gz_tbi, failOnDuplicate:true, failOnMismatch:true)
        .mix(ch_all_beds.no_bed)
        .set { ch_beds }

    //
    // Call the variants
    //

    BAM_STRUCTURAL_VARIANT_CALLING(
        ch_inputs.crams,
        ch_beds,
        ch_fasta,
        ch_fai,
        ch_bwa_index
    )

    ch_versions = ch_versions.mix(BAM_STRUCTURAL_VARIANT_CALLING.out.versions)
    ch_reports  = ch_reports.mix(BAM_STRUCTURAL_VARIANT_CALLING.out.reports)

    //
    // Genotype the variants
    //

    VCF_GENOTYPE_SV_PARAGRAPH(
        BAM_STRUCTURAL_VARIANT_CALLING.out.vcfs,
        ch_inputs.crams,
        ch_fasta,
        ch_fai
    )
    ch_versions = ch_versions.mix(VCF_GENOTYPE_SV_PARAGRAPH.out.versions)

    //
    // Annotate using Ensembl VEP
    //

    if(params.annotate) {
        VCF_ANNOTATE_VEP_ANNOTSV_VCFANNO(
            VCF_GENOTYPE_SV_PARAGRAPH.out.genotyped_vcfs,
            ch_small_variants_ready,
            ch_fasta,
            ch_fai,
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
