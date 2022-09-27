//
// Ploidy
//

include { BUILD_PLOIDY_MATRIX   } from '../../../modules/local/ploidy/build-ploidy-matrix/main'
include { PLOIDY_SCORE          } from '../../../modules/local/ploidy/ploidy-score/main'

workflow PLOIDY {
    take:
        bincov_matrix             // channel: [mandatory] [ bincov_matrix ] => The bin coverage matrix

    main:

    ch_versions = Channel.empty()

    BUILD_PLOIDY_MATRIX(
        bincov_matrix
    )

    ch_versions     = ch_versions.mix(BUILD_PLOIDY_MATRIX.out.versions)

    PLOIDY_SCORE(
        ploidy_matrix
    )

    ch_versions     = ch_versions.mix(PLOIDY_SCORE.out.versions)



    emit:
    ploidy_matrix   = BUILD_PLOIDY_MATRIX.out.ploidy_matrix
    ploidy_plots    = PLOIDY_SCORE.out.ploidy_plots

    versions        = ch_versions
}