//
// Gather Batch Evidence
//

// Import subworkflows
include { MAKE_BINCOV_MATRIX                            } from '../common-workflows/make-bincov-matrix'

// Import modules
include { GATK4_COLLECTREADCOUNTS as COLLECTREADCOUNTS  } from '../../../modules/nf-core/modules/gatk4/collectreadcounts/main'

workflow GATHER_BATCH_EVIDENCE {
    take:
        count_files             // channel: [mandatory] [ meta, count_file ] => The read counts of all samples
        bincov_matrix           // channel: [optional]  [ bincov_matrix ] => The bin coverage matrix
        bincov_matrix_index     // channel: [optional]  [ bincov_matrix_index ] => The bin coverage matrix index

    main:

    ch_versions = Channel.empty()

    if(!bincov_matrix && !bincov_matrix_index){
        MAKE_BINCOV_MATRIX(
            count_files
        )

        ch_versions = ch_versions.mix(MAKE_BINCOV_MATRIX.out.versions)

        bincov_matrix       = MAKE_BINCOV_MATRIX.out.merged_bincov
        bincov_matrix_index = MAKE_BINCOV_MATRIX.out.merged_bincov_index
    }

    ch_versions = Channel.empty()

    emit:
    versions            = ch_versions
}