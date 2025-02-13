//
// Run Delly
//

include { DELLY_CALL        } from '../../../modules/nf-core/delly/call/main'
include { BCFTOOLS_CONCAT   } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_SORT     } from '../../../modules/nf-core/bcftools/sort/main'
include { SVYNC             } from '../../../modules/nf-core/svync/main'

workflow BAM_VARIANT_CALLING_DELLY {
    take:
        ch_crams            // channel: [mandatory] [ meta, cram, crai ] => The aligned CRAMs per sample with the regions they should be called on
        ch_fasta            // channel: [mandatory] [ meta, fasta ] => The fasta reference file
        ch_fai              // channel: [mandatory] [ meta, fai ] => The index of the fasta reference file
        ch_svync_configs    // channel: [mandatory] [ configs ] => A list of svync config files

    main:

    def ch_versions     = Channel.empty()

    //
    // Calling variants using Delly
    //

    def sv_types = ["DEL", "INS", "INV", "DUP", "BND"]

    def ch_delly_input = ch_crams
        .combine(sv_types)
        .map { meta, cram, crai, sv_type ->
            def new_meta = meta + [caller:'delly', sv_type:sv_type]
            [ new_meta, cram, crai, [], [], [] ]
        }
        .dump(tag: 'delly_input', pretty: true)

    DELLY_CALL(
        ch_delly_input,
        ch_fasta,
        ch_fai
    )
    ch_versions = ch_versions.mix(DELLY_CALL.out.versions.first())

    def ch_concat_input = DELLY_CALL.out.bcf
        .join(DELLY_CALL.out.csi, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, bcf, csi ->
            def new_meta = meta - meta.subMap("sv_type")
            [ new_meta, bcf, csi ]
        }
        .groupTuple(size:5)

    BCFTOOLS_CONCAT(
        ch_concat_input
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions.first())

    BCFTOOLS_SORT(
        BCFTOOLS_CONCAT.out.vcf
    )
    ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions.first())

    def ch_delly_svync_config = ch_svync_configs
        .map { configs ->
            configs.find { config -> config.toString().contains("delly") }
        }

    def ch_delly_vcfs = BCFTOOLS_SORT.out.vcf
        .join(BCFTOOLS_SORT.out.tbi, failOnDuplicate:true, failOnMismatch:true)

    def ch_svync_input = ch_delly_vcfs
        .combine(ch_delly_svync_config)
        .dump(tag: 'delly_vcfs', pretty: true)

    SVYNC(
        ch_svync_input
    )
    ch_versions = ch_versions.mix(SVYNC.out.versions.first())

    def ch_out_vcfs = SVYNC.out.vcf
        .join(SVYNC.out.tbi, failOnDuplicate:true, failOnMismatch:true)

    emit:
    raw_vcfs    = ch_delly_vcfs // channel: [ val(meta), path(vcf), path(tbi) ]
    delly_vcfs  = ch_out_vcfs   // channel: [ val(meta), path(vcf), path(tbi) ]

    versions    = ch_versions
}
