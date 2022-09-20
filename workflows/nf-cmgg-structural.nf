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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { TABIX_BGZIPTABIX                  } from '../modules/nf-core/modules/tabix/bgziptabix/main'
include { BEDTOOLS_SORT                     } from '../modules/nf-core/modules/bedtools/sort/main'
include { GATK4_CREATESEQUENCEDICTIONARY    } from '../modules/nf-core/modules/gatk4/createsequencedictionary/main'
include { SAMTOOLS_FAIDX                    } from '../modules/nf-core/modules/samtools/faidx/main'
include { MULTIQC                           } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS       } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

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

    inputs = parse_input(ch_input).multiMap({ meta, cram, crai, bed=[], baseline=[] ->
        bed: [ meta, bed ]
        crams: [ meta, cram, crai ]
    })

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

    EVIDENCE_QC(
        GATHER_SAMPLE_EVIDENCE.out.vcfs,
        GATHER_SAMPLE_EVIDENCE.out.coverage_counts
    )

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
    Functions used in the pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def parse_input(input_csv) {

    // The samplesheet schema (change this to adjust the input check)
    def samplesheet_schema = [
        'columns': [
            'sample': [
                'content': 'meta',
                'meta_name': 'id,sample',
                'pattern': '',
            ],
            'cram': [
                'content': 'file',
                'pattern': '^.*\\.cram$',
            ],
            'crai': [
                'content': 'file',
                'pattern': '^.*\\.crai$',
            ],
            'bed': [
                'content': 'file',
                'pattern': '^.*\\.bed$',
            ],
        ],
        'required': ['sample','cram','crai'],
    ]    

    // Don't change these variables    
    def row_count = 1
    def all_columns = samplesheet_schema.columns.keySet().collect()
    def mandatory_columns = samplesheet_schema.required

    // Header checks
    Channel.value(input_csv).splitCsv(strip:true).first().map({ row ->

        if(row != all_columns) {
            def commons = all_columns.intersect(row)
            def diffs = all_columns.plus(row)
            diffs.removeAll(commons)

            if(diffs.size() > 0){
                def missing_columns = []
                def wrong_columns = []
                for(diff : diffs){
                    diff in all_columns ? missing_columns.add(diff) : wrong_columns.add(diff)
                }
                if(missing_columns.size() > 0){
                    exit 1, "[Samplesheet Error] The column(s) $missing_columns is/are not present. The header should look like: $all_columns"
                }
                else {
                    exit 1, "[Samplesheet Error] The column(s) $wrong_columns should not be in the header. The header should look like: $all_columns"
                }
            }
            else {
                exit 1, "[Samplesheet Error] The columns $row are not in the right order. The header should look like: $all_columns"
            }

        }
    })

    // Field checks + returning the channels
    Channel.value(input_csv).splitCsv(header:true, strip:true).map({ row ->

        row_count++

        // Check the mandatory columns
        def missing_mandatory_columns = []
        for(column : mandatory_columns) {
            row[column] ?: missing_mandatory_columns.add(column)
        }
        if(missing_mandatory_columns.size > 0){
            exit 1, "[Samplesheet Error] The mandatory column(s) $missing_mandatory_columns is/are empty on line $row_count"
        }

        def output = []
        def meta = [:]
        for(col : samplesheet_schema.columns) {
            key = col.key
            content = row[key]
            
            if(!(content ==~ col.value['pattern']) && col.value['pattern'] != '' && content != '') {
                exit 1, "[Samplesheet Error] The content of column '$key' on line $row_count does not match the pattern '${col.value['pattern']}'"
            }

            if(col.value['content'] == 'file'){
                output.add(content ? file(content, checkIfExists:true) : [])
            }
            else if(col.value['content'] == 'meta' && content != ''){
                for(meta_name : col.value['meta_name'].split(",")){
                    meta[meta_name] = content.replace(' ', '_')
                }
            }
        }

        output.add(0, meta)
        return output
    })
    
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
