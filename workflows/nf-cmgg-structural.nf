/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowNfcmggstructural.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
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

// Parse parameters
fasta           = params.fasta
fasta_fai       = params.fasta_fai
dict            = params.dict
allele_loci_vcf = params.allele_loci_vcf ?: []

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { GATHER_SAMPLE_EVIDENCE  } from '../subworkflows/local/gather-sample-evidence/main'
include { EVIDENCE_QC             } from '../subworkflows/local/evidence-QC/main'
include { GATHER_BATCH_EVIDENCE   } from '../subworkflows/local/gather-batch-evidence/main'

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
include { MULTIQC                           } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS       } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow NF_CMGG_STRUCTURAL {

    ch_versions = Channel.empty()
    ch_reports  = Channel.empty()

    //
    // Create optional inputs
    //

    if(!fasta_fai){
        SAMTOOLS_FAIDX(
            [ [], fasta ]
        )

        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
        fasta_fai   = SAMTOOLS_FAIDX.out.fai
    }

    if(!dict) {
        GATK4_CREATESEQUENCEDICTIONARY(
            fasta
        )

        ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
        dict        = GATK4_CREATESEQUENCEDICTIONARY.out.dict
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
        "bed"
    )

    ch_versions = ch_versions.mix(BEDTOOLS_SORT.out.versions)

    TABIX_BGZIPTABIX(
        BEDTOOLS_SORT.out.sorted
    )
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions)

    beds = BEDTOOLS_SORT.out.sorted.combine(TABIX_BGZIPTABIX.out.gz_tbi, by:0)

    //
    // Gather sample evidence
    //

    GATHER_SAMPLE_EVIDENCE(
        inputs.crams,
        beds,
        allele_loci_vcf,
        fasta,
        fasta_fai,
        dict
    )

    ch_versions = ch_versions.mix(GATHER_SAMPLE_EVIDENCE.out.versions)
    ch_reports  = ch_reports.mix(GATHER_SAMPLE_EVIDENCE.out.reports)

    //
    // Evidence QC
    //

    // EVIDENCE_QC(
    //     GATHER_SAMPLE_EVIDENCE.out.vcfs,
    //     GATHER_SAMPLE_EVIDENCE.out.coverage_counts,
    //     []
    // )

    // ch_versions = ch_versions.mix(EVIDENCE_QC.out.versions)

    //
    // Gather batch evidence
    //

    // GATHER_BATCH_EVIDENCE(
    //     GATHER_SAMPLE_EVIDENCE.out.coverage_counts,
    //     [], //EVIDENCE_QC.out.bincov_matrix,
    //     [], //EVIDENCE_QC.out.bincov_matrix_index
    //     [], // BAF files
    //     GATHER_SAMPLE_EVIDENCE.out.read_pairs,
    //     GATHER_SAMPLE_EVIDENCE.out.split_reads,
    //     GATHER_SAMPLE_EVIDENCE.out.site_depths,
    //     fasta,
    //     fasta_fai,
    //     dict
    // )

    // ch_versions = ch_versions.mix(GATHER_BATCH_EVIDENCE.out.versions)

    //
    // Cluster Batch
    //

    // CLUSTER_BATCH()

    //
    // Generate Batch Metrics
    //

    // GENERATE_BATCH_METRICS()

    //
    // Filter Batch
    //

    // FILTER_BATCH()

    //
    // Merge Batch Sites
    //

    // MERGE_BATCH_SITES()

    //
    // Genotype Batch
    //

    // GENOTYPE_BATCH()

    //
    // Make Cohort VCF
    //

    // MAKE_COHORT_VCF()

    //
    // Dump the software versions
    //

    CUSTOM_DUMPSOFTWAREVERSIONS(
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    ch_versions_yaml = CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect()

    //
    // Perform multiQC on all QC data
    //

    ch_multiqc_files = Channel.empty()

    ch_multiqc_files = ch_multiqc_files.mix(
                                        ch_versions_yaml,
                                        ch_reports.collect(),
                                        ch_multiqc_custom_config
                                        )

    MULTIQC (
        ch_multiqc_files.collect()
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
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
