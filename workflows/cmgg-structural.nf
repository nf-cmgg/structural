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
    params.fasta_fai,
    params.dict,
    params.allele_loci_vcf
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input, checkIfExists: true) } else { exit 1, 'Input samplesheet not specified!' }

// Check callers
def availableCallers = [
    "delly",
    "whamg",
    "manta",
    "gridss",
    "smoove"
]

for (caller in params.callers.tokenize(",")) {
    if(!(caller in availableCallers)) { exit 1, "The caller '${caller}' is not supported please specify a comma delimited list with on or more of the following callers: ${availableCallers}".toString() }
}

// Parse parameters
fasta           = Channel.fromPath(params.fasta).collect()
fasta_fai       = params.fasta_fai ? Channel.fromPath(params.fasta_fai).collect() : null
dict            = params.dict ? Channel.fromPath(params.dict).collect() : null
bwa_index       = params.bwa ? Channel.fromPath(params.bwa).map {[[],it]}.collect() : null
allele_loci_vcf = params.allele_loci_vcf ? Channel.fromPath(params.allele_loci_vcf).collect() : []

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
include { VCF_GENOTYPE_SV_DELLY             } from '../subworkflows/local/vcf_genotype_sv_delly/main'

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
include { GATK4_CREATESEQUENCEDICTIONARY    } from '../modules/nf-core/gatk4/createsequencedictionary/main'
include { SAMTOOLS_FAIDX                    } from '../modules/nf-core/samtools/faidx/main'
include { BWA_INDEX                         } from '../modules/nf-core/bwa/index/main'
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
    // Create optional inputs
    //

    if(!fasta_fai){
        SAMTOOLS_FAIDX(
            fasta.map {[[],it]}
        )

        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
        fasta_fai   = SAMTOOLS_FAIDX.out.fai.map { it[1] }.collect()
    }

    if(!dict) {
        GATK4_CREATESEQUENCEDICTIONARY(
            fasta
        )

        ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
        dict        = GATK4_CREATESEQUENCEDICTIONARY.out.dict.collect()
    }

    if(!bwa_index && params.callers.contains("gridss")){
        BWA_INDEX(
            fasta.map {[[id:'bwa'],it]}
        )

        ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
        bwa_index = BWA_INDEX.out.index.collect()
    }

    //
    // Create the input channel
    //

    SamplesheetConversion.convert(ch_input, file("${projectDir}/assets/schema_input.json"))
        .multiMap({ meta, cram, crai, bed, oed ->
            bed: [ meta, bed ]
            crams: [ meta, cram, crai ]
        })
        .set { inputs }

    //
    // Prepare the BED files
    //

    BEDTOOLS_SORT(
        inputs.bed,
        []
    )

    ch_versions = ch_versions.mix(BEDTOOLS_SORT.out.versions)

    TABIX_BGZIPTABIX(
        BEDTOOLS_SORT.out.sorted
    )
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions)

    beds = BEDTOOLS_SORT.out.sorted.combine(TABIX_BGZIPTABIX.out.gz_tbi, by:0)

    //
    // Call the variants
    //

    BAM_STRUCTURAL_VARIANT_CALLING(
        inputs.crams,
        beds,
        allele_loci_vcf,
        fasta,
        fasta_fai,
        dict,
        bwa_index
    )

    ch_versions = ch_versions.mix(BAM_STRUCTURAL_VARIANT_CALLING.out.versions)
    ch_reports  = ch_reports.mix(BAM_STRUCTURAL_VARIANT_CALLING.out.reports)

    //
    // Genotype the variants
    //

    VCF_GENOTYPE_SV_DELLY(
        BAM_STRUCTURAL_VARIANT_CALLING.out.vcfs,
        inputs.crams,
        fasta,
        fasta_fai
    )

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
