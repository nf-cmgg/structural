//
// Run Delly
//

include { DELLY_CALL        } from '../../../modules/nf-core/delly/call/main'
include { TABIX_TABIX       } from '../../../modules/nf-core/tabix/tabix/main'
include { SVYNC             } from '../../../modules/nf-core/svync/main'

workflow BAM_VARIANT_CALLING_DELLY {
    take:
        ch_crams            // channel: [mandatory] [ meta, cram, crai ] => The aligned CRAMs per sample with the regions they should be called on
        ch_fasta            // channel: [mandatory] [ meta, fasta ] => The fasta reference file
        ch_fai              // channel: [mandatory] [ meta, fai ] => The index of the fasta reference file
        ch_svync_configs    // channel: [mandatory] [ configs ] => A list of svync config files

    main:

    ch_versions     = Channel.empty()

    //
    // Calling variants using Delly
    //

    def ch_delly_input = ch_crams
        .map { meta, cram, crai ->
            [ meta, cram, crai, [], [], [] ]
        }
        .dump(tag: 'delly_input', pretty: true)

    DELLY_CALL(
        ch_delly_input,
        ch_fasta,
        ch_fai
    )
    ch_versions = ch_versions.mix(DELLY_CALL.out.versions.first())

    def ch_delly_svync_config = ch_svync_configs
        .map { configs ->
            configs.find { config -> config.toString().contains("delly") }
        }

    def ch_delly_vcfs = DELLY_CALL.out.bcf
        .join(DELLY_CALL.out.csi, failOnDuplicate:true, failOnMismatch:true)
        .combine(ch_delly_svync_config)
        .map { meta, vcf, tbi, config ->
            def new_meta = meta + [caller:"delly"]
            [ new_meta, vcf, tbi, config ]
        }
        .dump(tag: 'delly_vcfs', pretty: true)

    SVYNC(
        ch_delly_vcfs
    )
    ch_versions = ch_versions.mix(SVYNC.out.versions.first())

    TABIX_TABIX(
        SVYNC.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    def ch_out_vcfs = SVYNC.out.vcf
        .join(TABIX_TABIX.out.tbi)

    emit:
    delly_vcfs  = ch_out_vcfs // channel: [ val(meta), path(vcf), path(tbi) ]

    versions    = ch_versions
}
