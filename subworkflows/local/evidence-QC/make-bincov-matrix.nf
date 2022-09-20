//
// Gather Sample Evidence
//
include { SET_BINS                      } from '../../../modules/local/bincov-matrix/set-bins/main'
include { MAKE_BINCOV_MATRIX_COLUMNS    } from '../../../modules/local/bincov-matrix/make-bincov-matrix-columns/main'

workflow MAKE_BINCOV_MATRIX {
    take:
        count_files             // channel: [mandatory] [ meta, tsv ] => The read counts collected in GatherSampleEvidence

    main:

    ch_versions = Channel.empty()

    SET_BINS(
        count_files.first()
    )

    ch_versions = ch_versions.mix(SET_BINS.out.versions)

    MAKE_BINCOV_MATRIX_COLUMNS(
        count_files,
        SET_BINS.out.binsize,
        SET_BINS.out.bin_locs
    )

    ch_versions = ch_versions.mix(MAKE_BINCOV_MATRIX_COLUMNS.out.versions)
    MAKE_BINCOV_MATRIX_COLUMNS.out.bincov.view()

    emit:
    versions = ch_versions
}