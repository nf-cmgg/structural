//
// Estimate repeat size with Expanionhunter
//

include { EXPANSIONHUNTER } from '../../../modules/nf-core/expansionhunter/main'
include { BCFTOOLS_FILTER } from '../../../modules/nf-core/bcftools/filter/main'
include { TABIX_TABIX     } from '../../../modules/nf-core/tabix/tabix/main'

workflow BAM_REPEAT_ESTIMATION_EXPANSIONHUNTER {
    take:
        ch_crams        // channel: [mandatory] [ meta, cram, crai, bed ] => The aligned CRAMs per sample with the regions they should be called on
        ch_fasta        // channel: [mandatory] [ meta, fasta ] => The fasta reference file
        ch_fai          // channel: [mandatory] [ meta, fai ] => The index of the fasta reference file
        ch_catalogue    // channel: [mandatory] [ meta, catalogue ] => The expansionhunter catalogue

    main:

    ch_versions = Channel.empty()

    EXPANSIONHUNTER(
        ch_crams,
        ch_fasta,
        ch_fai,
        ch_catalogue
    )
    ch_versions = ch_versions.mix(EXPANSIONHUNTER.out.versions.first())

    BCFTOOLS_FILTER(
        EXPANSIONHUNTER.out.vcf
    )
    ch_versions = ch_versions.mix(BCFTOOLS_FILTER.out.versions.first())

    TABIX_TABIX(
        BCFTOOLS_FILTER.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    BCFTOOLS_FILTER.out.vcf
        .join(TABIX_TABIX.out.tbi, failOnDuplicate:true, failOnMismatch:true)
        .set { ch_vcfs }

    emit:
    vcfs                = ch_vcfs    // channel: [ val(meta), path(vcf), path(tbi) ]

    versions            = ch_versions
}
