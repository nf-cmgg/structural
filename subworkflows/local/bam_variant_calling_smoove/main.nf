//
// Run Delly
//

include { SMOOVE_CALL        } from '../../../modules/nf-core/smoove/call/main'
include { BCFTOOLS_SORT      } from '../../../modules/nf-core/bcftools/sort/main'
include { TABIX_TABIX        } from '../../../modules/nf-core/tabix/tabix/main'
include { SVYNC              } from '../../../modules/nf-core/svync/main'

workflow BAM_VARIANT_CALLING_SMOOVE {
    take:
        ch_crams    // channel: [mandatory] [ meta, cram, crai ] => The aligned CRAMs per sample with the regions they should be called on
        ch_fasta    // channel: [mandatory] [ meta, fasta ] => The fasta reference file
        ch_fai      // channel: [mandatory] [ meta, fai ] => The index of the fasta reference file

    main:

    ch_versions     = Channel.empty()

    //
    // Calling variants using Smoove
    //

    ch_crams
        .map { meta, cram, crai ->
            [ meta, cram, crai, [] ]
        }
        .dump(tag: 'smoove_input', pretty: true)
        .set { ch_smoove_input }

    SMOOVE_CALL(
        ch_smoove_input,
        ch_fasta.map{it[1]},
        ch_fai.map{it[1]}
    )

    ch_versions = ch_versions.mix(SMOOVE_CALL.out.versions.first())

    BCFTOOLS_SORT(
        SMOOVE_CALL.out.vcf
    )
    ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions.first())

    BCFTOOLS_SORT.out.vcf
        .map(
            { meta, vcf ->
                new_meta = meta + [caller:'smoove']
                [ new_meta, vcf, [] ]
            }
        )
        .dump(tag: 'smoove_vcfs', pretty: true)
        .set { ch_smoove_vcfs }

    Channel.fromPath("${projectDir}/assets/svync/smoove.yaml")
        .map { [[], it] }
        .collect()
        .set { ch_svync_config }

    SVYNC(
        ch_smoove_vcfs,
        ch_svync_config
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
    smoove_vcfs = ch_out_vcfs    // channel: [ val(meta), path(vcf), path(tbi) ]

    versions    = ch_versions
}
