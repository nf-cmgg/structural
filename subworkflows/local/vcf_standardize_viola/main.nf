//
// Standardize with viola
//

include { VIOLA                                         } from '../../../modules/local/viola/main'

workflow VCF_STANDARDIZE_VIOLA {
    take:
        ch_vcfs        // channel: [mandatory] [ meta, vcf ] => The called VCFs containing SVs

    main:

    ch_versions     = Channel.empty()

    VIOLA(
        ch_vcfs
    )
    ch_versions = ch_versions.mix(VIOLA.out.versions.first())

    VIOLA.out.vcf.set { ch_standardized_vcfs}

    emit:
    standardized_vcfs = ch_standardized_vcfs // channel: [ val(meta), path(vcf) ]
    versions          = ch_versions
}
