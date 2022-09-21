//
// Evidence QC
//

// Import subworkflows
include { MAKE_BINCOV_MATRIX    } from './make-bincov-matrix'

// Import modules
include { CALCMEDCOV            } from '../../../modules/local/calcmedcov/main'
include { BUILD_PLOIDY_MATRIX   } from '../../../modules/local/ploidy/build-ploidy-matrix/main'
include { PLOIDY_SCORE          } from '../../../modules/local/ploidy/ploidy-score/main'
include { WGD_SCORE             } from '../../../modules/local/WGD-score/main'

include { BEDTOOLS_INTERSECT    } from '../../../modules/nf-core/modules/bedtools/intersect/main'

workflow EVIDENCE_QC {
    take:
        vcfs                    // channel: [optional]  [ meta, vcf, tbi ] => The called VCFs from GatherSampleEvidence
        count_files             // channel: [mandatory] [ meta, tsv ] => The read counts collected in GatherSampleEvidence
        wgd_scoring_mask        // channel: [optional]  [ wgd_scoring_mask ] => The WGD scoring mask BED file

    main:

    ch_versions = Channel.empty()

    MAKE_BINCOV_MATRIX(
        count_files
    )

    ch_versions = ch_versions.mix(MAKE_BINCOV_MATRIX.out.versions)

    bincov_matrix       = MAKE_BINCOV_MATRIX.out.merged_bincov
    bincov_matrix_gz    = MAKE_BINCOV_MATRIX.out.merged_bincov_gz

    // FIX THIS!!
    // CALCMEDCOV(
    //     bincov_matrix
    // ).median_cov_file.view()


    //
    // Ploidy
    //

    if(params.run_ploidy){
            
        BUILD_PLOIDY_MATRIX(
            bincov_matrix_gz
        )

        ch_versions.mix(BUILD_PLOIDY_MATRIX.out.versions)

        PLOIDY_SCORE(
            BUILD_PLOIDY_MATRIX.out.ploidy_matrix
        )

        ch_versions = ch_versions.mix(PLOIDY_SCORE.out.versions)
    
    }

    //
    // WGD
    //

    if(wgd_scoring_mask){

        BEDTOOLS_INTERSECT(
            bincov_matrix.combine(wgd_scoring_mask).map({ bincov, wgd_mask -> [ [], bincov, wgd_mask ]}),
            'bed'
        )

        ch_versions = ch_versions.mix(BEDTOOLS_INTERSECT.out.versions)

        WGD_SCORE(
            BEDTOOLS_INTERSECT.out.intersect.combine(wgd_scoring_mask)
        )

    }

    emit:
    versions            = ch_versions
}