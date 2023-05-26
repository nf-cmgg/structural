//
// Run Whamg
//
include { WHAMG                       } from '../../../modules/nf-core/whamg/main'
include { SAMTOOLS_CONVERT            } from '../../../modules/nf-core/samtools/convert/main'
include { TABIX_TABIX as TABIX_WHAMG  } from '../../../modules/nf-core/tabix/tabix/main'

workflow BAM_VARIANT_CALLING_WHAMG {
    take:
        crams                   // channel: [mandatory] [ meta, cram, crai ] => The aligned CRAMs per sample with the regions they should be called on
        fasta                   // channel: [mandatory] [ fasta ] => The fasta reference file
        fai                     // channel: [mandatory] [ fai ] => The index of the fasta reference file

    main:

    ch_versions      = Channel.empty()

    //
    // Convert the CRAMs to BAMs
    //

    SAMTOOLS_CONVERT(
        crams,
        fasta,
        fai
    )

    SAMTOOLS_CONVERT.out.alignment_index.set { bams }
    ch_versions = ch_versions.mix(SAMTOOLS_CONVERT.out.versions)

    //
    // Calling variants using Whamg
    //

    WHAMG(
        bams,
        fasta,
        fai
    )

    ch_versions = ch_versions.mix(WHAMG.out.versions)

    WHAMG.out.vcf
        .join(WHAMG.out.tbi, failOnMismatch:true, failOnDuplicate:true)
        .set { whamg_vcfs }

    whamg_vcfs
        .map(
            { meta, vcf, tbi ->
                new_meta = meta.findAll {true} + [caller:"whamg"]
                [ new_meta, vcf, tbi ]
            }
        )
        .dump(tag: 'whamg_vcfs', pretty: true)
        .set { whamg_vcfs }

    emit:
    whamg_vcfs
    versions = ch_versions
}
