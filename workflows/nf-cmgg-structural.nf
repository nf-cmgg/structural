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
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input, checkIfExists: true) } else { exit 1, 'Input samplesheet not specified!' }

// Parse parameters
fasta     = params.fasta
fasta_fai = params.fasta_fai
dict      = params.dict

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
include { GATHER_SAMPLE_EVIDENCE      } from '../subworkflows/local/gather-sample-evidence'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { TABIX_BGZIPTABIX            } from '../modules/nf-core/modules/tabix/bgziptabix/main'
include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow NF_CMGG_STRUCTURAL {

    ch_versions = Channel.empty()

    inputs = parse_input(ch_input).multiMap({ meta, cram, crai, bed ->
        bed: [ meta, bed ]
        crams: [ meta, cram, crai ]
    })

    TABIX_BGZIPTABIX(
        inputs.bed
    )

    beds = inputs.bed.combine(TABIX_BGZIPTABIX.out.gz_tbi, by:0)

    GATHER_SAMPLE_EVIDENCE(
        inputs.crams,
        beds,
        fasta,
        fasta_fai,
        dict
    )

    GATHER_SAMPLE_EVIDENCE.out.vcfs.view()


    // MULTIQC (
    //     ch_multiqc_files.collect()
    // )
    // multiqc_report = MULTIQC.out.report.toList()
    // ch_versions    = ch_versions.mix(MULTIQC.out.versions)
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
                'meta_name': 'id',
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
            ]
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
            
            if(!(content ==~ col.value['pattern']) && col.value['pattern'] != '') {
                exit 1, "[Samplesheet Error] The content of column '$key' on line $row_count does not match the pattern '${col.value['pattern']}'"
            }

            if(col.value['content'] == 'file'){
                output.add(content ? file(content, checkIfExists:true) : [])
            }
            else if(col.value['content'] == 'meta'){
                meta[col.value['meta_name']] = content
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
