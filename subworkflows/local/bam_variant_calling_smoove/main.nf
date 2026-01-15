//
// Run Smoove
//

include { SMOOVE_CALL                   } from '../../../modules/nf-core/smoove/call/main'
include { BCFTOOLS_SORT                 } from '../../../modules/nf-core/bcftools/sort/main'
include { TABIX_TABIX                   } from '../../../modules/nf-core/tabix/tabix/main'
include { SVYNC                         } from '../../../modules/nf-core/svync/main'

workflow BAM_VARIANT_CALLING_SMOOVE {
    take:
        ch_crams            // channel: [mandatory] [ meta, cram, crai ] => The aligned CRAMs per sample with the regions they should be called on
        ch_fasta            // channel: [mandatory] [ meta, fasta ] => The fasta reference file
        ch_fai              // channel: [mandatory] [ meta, fai ] => The index of the fasta reference file
        ch_svync_configs    // channel: [mandatory] [ configs ] => A list of svync config files

    main:

    def ch_versions     = channel.empty()

    //
    // Calling variants using Smoove
    //

    def ch_smoove_input = ch_crams
        .map { meta, cram, crai ->
            def new_meta = meta + [caller:'smoove']
            [ new_meta, cram, crai, [] ]
        }
        .dump(tag: 'smoove_input', pretty: true)

    SMOOVE_CALL(
        ch_smoove_input,
        ch_fasta,
        ch_fai
    )

    ch_versions = ch_versions.mix(SMOOVE_CALL.out.versions.first())

    BCFTOOLS_SORT(
        SMOOVE_CALL.out.vcf
    )
    ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions.first())

    def ch_smoove_svync_config = ch_svync_configs
        .map { configs ->
            configs.find { config -> config.toString().contains("smoove") }
        }

    def ch_smoove_vcfs = BCFTOOLS_SORT.out.vcf
        .join(BCFTOOLS_SORT.out.tbi)

    def ch_svync_input = ch_smoove_vcfs
        .combine(ch_smoove_svync_config)
        .dump(tag: 'smoove_vcfs', pretty: true)

    SVYNC(
        ch_svync_input
    )
    ch_versions = ch_versions.mix(SVYNC.out.versions.first())

    TABIX_TABIX(
        SVYNC.out.vcf
    )

    def ch_out_vcfs = SVYNC.out.vcf
        .join(TABIX_TABIX.out.index)

    emit:
    raw_vcfs    = ch_smoove_vcfs    // channel: [ val(meta), path(vcf), path(tbi) ]
    smoove_vcfs = ch_out_vcfs       // channel: [ val(meta), path(vcf), path(tbi) ]

    versions    = ch_versions
}
