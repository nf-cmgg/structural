//
// Annotate the VCFs
//

include { BCFTOOLS_CONCAT   } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_SORT     } from '../../../modules/nf-core/bcftools/sort/main'

workflow VCF_CONCAT_BCFTOOLS {
    take:
        ch_vcfs         // channel: [mandatory] [ val(meta), path(vcf), path(tbi) ] VCFs containing the called structural variants
        val_count_types // value:   [mandatory] => The amount of different variant types

    main:

    def ch_versions = channel.empty()

    def ch_concat_input = ch_vcfs
        .groupTuple(size:val_count_types)

    BCFTOOLS_CONCAT(ch_concat_input)
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions.first())

    BCFTOOLS_SORT(BCFTOOLS_CONCAT.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions.first())

    def ch_concat_vcfs = BCFTOOLS_SORT.out.vcf
        .join(BCFTOOLS_SORT.out.tbi, failOnDuplicate:true, failOnMismatch:true)

    emit:
    vcfs            = ch_concat_vcfs  // channel: [ val(meta), path(vcf), path(tbi) ]

    versions        = ch_versions
}
