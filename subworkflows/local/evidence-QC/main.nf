//
// Evidence QC
//

// Import subworkflows
include { MAKE_BINCOV_MATRIX } from './make-bincov-matrix'

// Import modules

workflow EVIDENCE_QC {
    take:
        vcfs                    // channel: [optional]  [ meta, vcf, tbi ] => The called VCFs from GatherSampleEvidence
        count_files             // channel: [mandatory] [ meta, tsv ] => The read counts collected in GatherSampleEvidence

    main:

    ch_versions = Channel.empty()

    MAKE_BINCOV_MATRIX(
        count_files
    )


    emit:
    versions            = ch_versions
}