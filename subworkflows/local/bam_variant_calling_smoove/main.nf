//
// Run Delly
//

include { REVERSE_BED        } from '../../../modules/local/reversebed/main'

include { SMOOVE_CALL        } from '../../../modules/nf-core/smoove/call/main'
include { TABIX_TABIX        } from '../../../modules/nf-core/tabix/tabix/main'

workflow BAM_VARIANT_CALLING_SMOOVE {
    take:
        crams                   // channel: [mandatory] [ meta, cram, crai ] => The aligned CRAMs per sample with the regions they should be called on
        beds                    // channel: [optional]  [ meta, bed, bed_gz, bed_gz_tbi ] => A channel containing the normal BED, the bgzipped BED and its index file
        fasta                   // channel: [mandatory] [ fasta ] => The fasta reference file
        fasta_fai               // channel: [mandatory] [ fasta_fai ] => The index of the fasta reference file

    main:

    ch_versions     = Channel.empty()

    //
    // Reverse the BED file (It will only contain the regions that aren't of interest now)
    //

    beds
        .map(
            { meta, bed, bed_gz, bed_gz_tbi ->
                [ meta, bed ]
            }
        )
        .set { reverse_input }

    REVERSE_BED(
        reverse_input,
        fasta_fai
    )

    ch_versions = ch_versions.mix(REVERSE_BED.out.versions)

    //
    // Calling variants using Smoove
    //

    crams
        .join(REVERSE_BED.out.bed)
        .dump(tag: 'smoove_input', pretty: true)
        .set { smoove_input }

    SMOOVE_CALL(
        smoove_input,
        fasta,
        fasta_fai
    )

    ch_versions = ch_versions.mix(SMOOVE_CALL.out.versions)

    TABIX_TABIX(
        SMOOVE_CALL.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

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
        .set { smoove_vcfs }

    emit:
    smoove_vcfs
    versions = ch_versions
}
