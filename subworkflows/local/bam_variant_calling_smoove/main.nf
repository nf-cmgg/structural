//
// Run Delly
//

include { SMOOVE_CALL        } from '../../../modules/nf-core/smoove/call/main'
include { BCFTOOLS_SORT      } from '../../../modules/nf-core/bcftools/sort/main'
include { TABIX_TABIX        } from '../../../modules/nf-core/tabix/tabix/main'

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

    TABIX_TABIX(
        BCFTOOLS_SORT.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    BCFTOOLS_SORT.out.vcf
        .combine(TABIX_TABIX.out.tbi, by:0)
        .map(
            { meta, vcf, tbi ->
                new_meta = meta.clone()
                new_meta.caller = "smoove"
                [ new_meta, vcf, tbi ]
            }
        )
        .dump(tag: 'smoove_vcfs', pretty: true)
        .set { ch_smoove_vcfs }

    emit:
    smoove_vcfs = ch_smoove_vcfs    // channel: [ val(meta), path(vcf), path(tbi) ]

    versions    = ch_versions
}
