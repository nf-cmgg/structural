//
// Run Manta
//
include { MANTA_GERMLINE         } from '../../../modules/nf-core/manta/germline/main'
include { MANTA_CONVERTINVERSION } from '../../../modules/nf-core/manta/convertinversion/main'

workflow BAM_VARIANT_CALLING_MANTA {
    take:
        ch_crams    // channel: [mandatory] [ meta, cram, crai ] => The aligned CRAMs per sample with the regions they should be called on
        ch_fasta    // channel: [mandatory] [ meta, fasta ] => The fasta reference file
        ch_fai      // channel: [mandatory] [ meta, fai ] => The index of the fasta reference file

    main:

    ch_versions     = Channel.empty()

    //
    // Calling variants using Manta
    //

    ch_crams
        .map { meta, cram, crai ->
            [ meta, cram, crai, [], [] ]
        }
        .dump(tag: 'manta_input', pretty: true)
        .set { ch_manta_input }

    MANTA_GERMLINE(
        ch_manta_input,
        ch_fasta,
        ch_fai
    )

    ch_versions = ch_versions.mix(MANTA_GERMLINE.out.versions.first())

    //
    // Reformat the inversions into single inverted sequence junctions
    //

    MANTA_CONVERTINVERSION(
        MANTA_GERMLINE.out.diploid_sv_vcf,
        ch_fasta.map{it[1]}
    )

    ch_versions = ch_versions.mix(MANTA_CONVERTINVERSION.out.versions.first())

    MANTA_CONVERTINVERSION.out.vcf
        .join(MANTA_CONVERTINVERSION.out.tbi, failOnDuplicate:true, failOnMismatch:true)
        .map(
            { meta, vcf, tbi ->
                new_meta = meta.clone()
                new_meta.caller = "manta"
                [ new_meta, vcf, tbi ]
            }
        )
        .dump(tag: 'manta_vcfs', pretty: true)
        .set { ch_manta_vcfs }

    emit:
    manta_vcfs  = ch_manta_vcfs  // channel: [ val(meta), path(vcf), path(tbi) ]

    versions    = ch_versions
}
