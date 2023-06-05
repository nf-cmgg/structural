//
// Standardize with viola
//

include { VIOLA                                         } from '../../../modules/local/viola/main'
include { BCFTOOLS_FILTER                               } from '../../../modules/nf-core/bcftools/filter/main'

workflow VCF_STANDARDIZE_VIOLA {
    take:
        ch_vcfs        // channel: [mandatory] [ meta, vcf ] => The called VCFs containing SVs

    main:

    ch_versions     = Channel.empty()

    VIOLA(
        ch_vcfs
    )
    ch_versions = ch_versions.mix(VIOLA.out.versions.first())

    VIOLA.out.vcf
        .branch { meta, vcf ->
            gridss: meta.caller == "gridss"
            no_gridss: meta.caller != "gridss"
        }
        .set { ch_viola_output }

    BCFTOOLS_FILTER(
        ch_viola_output.gridss
    )
    ch_versions = ch_versions.mix(BCFTOOLS_FILTER.out.versions.first())

    ch_viola_output.no_gridss
        .mix(BCFTOOLS_FILTER.out.vcf)
        .set { ch_standardized_vcfs}

    emit:
    standardized_vcfs = ch_standardized_vcfs // channel: [ val(meta), path(vcf) ]
    versions          = ch_versions
}
