//
// Annotate the VCFs
//

include { BCFTOOLS_CONCAT   } from '../../../modules/nf-core/bcftools/concat/main'
include { TABIX_TABIX       } from '../../../modules/nf-core/tabix/tabix/main'

workflow VCF_CONCAT_BCFTOOLS {
    take:
        ch_vcfs         // channel: [mandatory] [ val(meta), path(vcf), path(tbi) ] VCFs containing the called structural variants
        val_count_types // value:   [mandatory] => The amount of different variant types

    main:

    ch_versions = Channel.empty()
    
    ch_vcfs
        .groupTuple(size:val_count_types)
        .set { ch_concat_input }

    BCFTOOLS_CONCAT(ch_concat_input)
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions.first())

    TABIX_TABIX(BCFTOOLS_CONCAT.out.vcf)
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    BCFTOOLS_CONCAT.out.vcf
        .join(TABIX_TABIX.out.tbi, failOnDuplicate:true, failOnMismatch:true)
        .set { ch_concat_vcfs }

    emit:
    vcfs            = ch_concat_vcfs  // channel: [ val(meta), path(vcf), path(tbi) ]

    versions        = ch_versions
}
