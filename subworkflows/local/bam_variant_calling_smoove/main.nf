//
// Run Delly
//

include { REVERSE_BED        } from '../../../modules/local/reversebed/main'

include { SMOOVE_CALL        } from '../../../modules/nf-core/smoove/call/main'
include { TABIX_TABIX        } from '../../../modules/nf-core/tabix/tabix/main'

workflow BAM_VARIANT_CALLING_SMOOVE {
    take:
        ch_crams    // channel: [mandatory] [ meta, cram, crai ] => The aligned CRAMs per sample with the regions they should be called on
        ch_beds     // channel: [optional]  [ meta, bed, bed_gz, bed_gz_tbi ] => A channel containing the normal BED, the bgzipped BED and its index file
        ch_fasta    // channel: [mandatory] [ fasta ] => The fasta reference file
        ch_fai      // channel: [mandatory] [ fai ] => The index of the fasta reference file

    main:

    ch_versions     = Channel.empty()

    //
    // Reverse the BED file (It will only contain the regions that aren't of interest now)
    //

    ch_beds
        .branch { meta, bed, bed_gz, bed_gz_tbi ->
            bed: bed
                return [ meta, bed ]
            no_bed: !bed
                return [ meta, [] ]
        }
        .set { ch_reverse_input }

    REVERSE_BED(
        ch_reverse_input.bed,
        ch_fai
    )

    ch_versions = ch_versions.mix(REVERSE_BED.out.versions.first())

    REVERSE_BED.out.bed
        .mix(ch_reverse_input.no_bed)
        .set { ch_reversed_beds }

    //
    // Calling variants using Smoove
    //

    ch_crams
        .join(ch_reversed_beds, failOnMismatch:true, failOnDuplicate:true)
        .dump(tag: 'smoove_input', pretty: true)
        .set { ch_smoove_input }

    SMOOVE_CALL(
        ch_smoove_input,
        ch_fasta,
        ch_fai
    )

    ch_versions = ch_versions.mix(SMOOVE_CALL.out.versions.first())

    TABIX_TABIX(
        SMOOVE_CALL.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    SMOOVE_CALL.out.vcf
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
