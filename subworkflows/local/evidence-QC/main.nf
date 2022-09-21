//
// Evidence QC
//

// Import subworkflows
include { MAKE_BINCOV_MATRIX    } from './make-bincov-matrix'

// Import modules
include { CALCMEDCOV            } from '../../../modules/local/calcmedcov/main'
include { BUILD_PLOIDY_MATRIX   } from '../../../modules/local/ploidy/build-ploidy-matrix/main'
include { PLOIDY_SCORE          } from '../../../modules/local/ploidy/ploidy-score/main'

workflow EVIDENCE_QC {
    take:
        vcfs                    // channel: [optional]  [ meta, vcf, tbi ] => The called VCFs from GatherSampleEvidence
        count_files             // channel: [mandatory] [ meta, tsv ] => The read counts collected in GatherSampleEvidence

    main:

    ch_versions = Channel.empty()

    MAKE_BINCOV_MATRIX(
        count_files
    )

    ch_versions = ch_versions.mix(MAKE_BINCOV_MATRIX.out.versions)

    bincov_matrix = MAKE_BINCOV_MATRIX.out.merged_bincov

    // FIX THIS!!
    // CALCMEDCOV(
    //     bincov_matrix
    // ).median_cov_file.view()

    if(params.run_ploidy){
            
        BUILD_PLOIDY_MATRIX(
            bincov_matrix
        )

        ch_versions.mix(BUILD_PLOIDY_MATRIX.out.versions)

        PLOIDY_SCORE(
            BUILD_PLOIDY_MATRIX.out.ploidy_matrix
        )

        ch_versions = ch_versions.mix(PLOIDY_SCORE.out.versions)
    
    }

    emit:
    versions            = ch_versions
}