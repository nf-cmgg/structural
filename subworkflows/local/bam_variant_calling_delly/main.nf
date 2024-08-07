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

    ch_crams
        .map { meta, cram, crai ->
            [ meta, cram, crai, [], [], [] ]
        }
        .dump(tag: 'delly_input', pretty: true)
        .set { ch_delly_input }

    DELLY_CALL(
        ch_delly_input,
        ch_fasta,
        ch_fai
    )
    ch_versions = ch_versions.mix(DELLY_CALL.out.versions.first())

    ch_svync_configs
        .map {
            it.find { it.toString().contains("delly") }
        }
        .set { ch_delly_svync_config }

    DELLY_CALL.out.bcf
        .join(DELLY_CALL.out.csi, failOnDuplicate:true, failOnMismatch:true)
        .combine(ch_delly_svync_config)
        .map { meta, vcf, tbi, config ->
            new_meta = meta + [caller:"delly"]
            [ new_meta, vcf, tbi, config ]
        }
        .dump(tag: 'delly_vcfs', pretty: true)
        .set { ch_delly_vcfs }

    SVYNC(
        ch_delly_vcfs
    )
    ch_versions = ch_versions.mix(SVYNC.out.versions.first())

    TABIX_TABIX(
        SVYNC.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    SVYNC.out.vcf
        .join(TABIX_TABIX.out.tbi)
        .set { ch_out_vcfs }

    emit:
    delly_vcfs  = ch_out_vcfs // channel: [ val(meta), path(vcf), path(tbi) ]

    versions    = ch_versions
}
