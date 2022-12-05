//
// Evidence QC
//

// Import subworkflows
include { MAKE_BINCOV_MATRIX    } from '../common-workflows/make-bincov-matrix'
include { PLOIDY                } from '../common-workflows/ploidy'

// Import modules
include { CALCMEDCOV            } from '../../../modules/local/calcmedcov/main'
include { WGD_SCORE             } from '../../../modules/local/WGD-score/main'
include { INDVIDUAL_QC          } from '../../../modules/local/vcf-QC/individual-QC/main'
include { PICK_OUTLIERS         } from '../../../modules/local/vcf-QC/pick-outliers/main'

include { BEDTOOLS_INTERSECT    } from '../../../modules/nf-core/bedtools/intersect/main'
include { TABIX_BGZIP           } from '../../../modules/nf-core/tabix/bgzip/main'

workflow EVIDENCE_QC {
    take:
        vcfs                    // channel: [optional]  [ meta, vcf, tbi ] => The called VCFs from GatherSampleEvidence
        count_files             // channel: [mandatory] [ meta, tsv ] => The read counts collected in GatherSampleEvidence
        wgd_scoring_mask        // channel: [optional]  [ wgd_scoring_mask ] => The WGD scoring mask BED file

    main:

    ch_versions         = Channel.empty()

    low_outliers        = Channel.empty()
    high_outliers       = Channel.empty()

    ploidy_matrix       = Channel.empty()
    ploidy_plots        = Channel.empty()

    wgd_dist            = Channel.empty()
    wgd_matrix          = Channel.empty()
    wgd_scores          = Channel.empty()

    bincov_median       = Channel.empty()


    MAKE_BINCOV_MATRIX(
        count_files
    )

    ch_versions = ch_versions.mix(MAKE_BINCOV_MATRIX.out.versions)

    MAKE_BINCOV_MATRIX.out.merged_bincov.set { bincov_matrix }
    MAKE_BINCOV_MATRIX.out.merged_bincov_index.set { bincov_matrix_index }

    // FIX THIS!!
    CALCMEDCOV(
        bincov_matrix
    ).median_cov_file.view()

    // bincov_median = bincov_median.mix(CALCMEDCOV.out.median_cov_file)


    //
    // Ploidy (BROKEN FOR NOW)
    //

    if(params.run_ploidy) {

        PLOIDY(
            bincov_matrix
        )

        ploidy_matrix   = ploidy_matrix.mix(PLOIDY.out.ploidy_matrix)
        ploidy_plots    = ploidy_plots.mix(PLOIDY.out.ploidy_plots)
        ch_versions     = ch_versions.mix(PLOIDY.out.versions)

    }

    //
    // WGD (untested, but should work)
    //

    if(wgd_scoring_mask) {

        BEDTOOLS_INTERSECT(
            bincov_matrix
                .combine(wgd_scoring_mask)
                .map(
                    { bincov, wgd_mask ->
                        [ [], bincov, wgd_mask ]
                    }
                ),
            'bed'
        )

        wgd_matrix  = wgd_matrix.mix(BEDTOOLS_INTERSECT.out.intersect)

        ch_versions = ch_versions.mix(BEDTOOLS_INTERSECT.out.versions)

        WGD_SCORE(
            wgd_matrix.combine(wgd_scoring_mask)
        )

        WGD_SCORE.out.dist.set { wgd_dist }
        WGD_SCORE.out.scores.set { wgd_scores }
        ch_versions = ch_versions.mix(WGD_SCORE.out.versions)

    }

    //
    // VCF QC (BROKEN FOR NOW)
    //

    if(params.run_vcf_qc) {

        INDVIDUAL_QC(
            vcfs.map({ meta, vcf, tbi -> [ meta, vcf ]})
        )

        ch_versions = ch_versions.mix(INDVIDUAL_QC.out.versions)

        INDVIDUAL_QC.out.stat
            .branch(
                { meta, stat ->
                    valid: stat.countLines() > 1
                    invalid: stat.countLines() == 1
                }
            )
            .set { pick_outliers_input }

        pick_outliers_input.valid.dump(tag: 'pick_outliers_input_valid', pretty: true)
        pick_outliers_input.invalid.dump(tag: 'pick_outliers_input_invalid', pretty: true)

        PICK_OUTLIERS(
            pick_outliers_input.valid
        )

        low_outliers    = low_outliers.mix(PICK_OUTLIERS.out.low)
        high_outliers   = high_outliers.mix(PICK_OUTLIERS.out.high)
        ch_versions     = ch_versions.mix(PICK_OUTLIERS.out.versions)

    }

    emit:
    high_outliers
    low_outliers

    ploidy_matrix
    ploidy_plots

    wgd_dist
    wgd_matrix
    wgd_scores

    bincov_matrix
    bincov_matrix_index
    bincov_median

    versions            = ch_versions
}
