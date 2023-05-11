//
// Run Expansionhunter
//

include { EXPANSIONHUNTER } from '../../../modules/nf-core/expansionhunter/main'
include { TABIX_TABIX     } from '../../../modules/nf-core/tabix/tabix/main'

workflow BAM_REPEATS_ESTIMATION_EXPANSIONHUNTER {
    take:
        ch_crams                // channel: [mandatory] [ meta, cram, crai ] => The aligned CRAMs per sample with the regions they should be called on
        ch_fasta                // channel: [mandatory] [ fasta ] => The fasta reference file
        ch_fai                  // channel: [mandatory] [ fai ] => The index of the fasta reference file
        ch_variant_catalogue    // channel: [mandatory] [ meta, variant_catalogue ] => The variant catalogue for expansionhunter

    main:

    ch_versions     = Channel.empty()

    EXPANSIONHUNTER(
        ch_crams,
        ch_fasta.map { [[], it] },
        ch_fai.map { [[], it] },
        ch_variant_catalogue
    )
    ch_versions = ch_versions.mix(EXPANSIONHUNTER.out.versions.first())

    TABIX_TABIX(
        EXPANSIONHUNTER.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    EXPANSIONHUNTER.out.vcf
        .join(TABIX_TABIX.out.tbi, failOnDuplicate:true, failOnMismatch:true)
        .set { ch_expansionhunter_vcfs }

    emit:
    expansionhunter_vcfs = ch_expansionhunter_vcfs // channel: [ val(meta), path(vcf), path(tbi) ]

    versions    = ch_versions
}
