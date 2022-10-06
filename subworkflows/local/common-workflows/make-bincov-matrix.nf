//
// Make Bincov Matrix
//

include { SET_BINS                      } from '../../../modules/local/bincov-matrix/set-bins/main'
include { MAKE_BINCOV_MATRIX_COLUMNS    } from '../../../modules/local/bincov-matrix/make-bincov-matrix-columns/main'
include { ZPASTE                        } from '../../../modules/local/bincov-matrix/zpaste/main'

workflow MAKE_BINCOV_MATRIX {
    take:
        count_files             // channel: [mandatory] [ meta, tsv ] => The read counts collected in GatherSampleEvidence

    main:

    ch_versions = Channel.empty()

    SET_BINS(
        count_files.first()
    )

    ch_versions = ch_versions.mix(SET_BINS.out.versions)

    make_bincov_matrix_input = count_files.combine(SET_BINS.out.binsize.splitText())
                                        .map({ meta, count_file, binsize ->
                                            new_meta = meta.clone()
                                            new_meta.binsize = binsize.replace("\n","")
                                            [ new_meta, count_file]
                                        })

    MAKE_BINCOV_MATRIX_COLUMNS(
        make_bincov_matrix_input,
        SET_BINS.out.bin_locs
    )

    ch_versions = ch_versions.mix(MAKE_BINCOV_MATRIX_COLUMNS.out.versions)

    zpaste_input = MAKE_BINCOV_MATRIX_COLUMNS.out.bincov
                    .map({ meta, bincov -> [ bincov ]})
                    .collect()

    ZPASTE(
        zpaste_input
    )

    ch_versions = ch_versions.mix(ZPASTE.out.versions)


    emit:
    merged_bincov          = ZPASTE.out.matrix_file
    merged_bincov_index    = ZPASTE.out.matrix_file_index

    versions               = ch_versions
}
