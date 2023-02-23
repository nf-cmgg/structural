//
// Run Manta
//
include { MANTA_GERMLINE         } from '../../../modules/nf-core/manta/germline/main'
include { MANTA_CONVERTINVERSION } from '../../../modules/nf-core/manta/convertinversion/main'
include { BCFTOOLS_REHEADER      } from '../../../modules/nf-core/bcftools/reheader/main'

workflow BAM_VARIANT_CALLING_MANTA {
    take:
        crams                   // channel: [mandatory] [ meta, cram, crai ] => The aligned CRAMs per sample with the regions they should be called on
        beds                    // channel: [optional]  [ meta, bed, bed_gz, bed_gz_tbi ] => A channel containing the normal BED, the bgzipped BED and its index file
        fasta                   // channel: [mandatory] [ fasta ] => The fasta reference file
        fasta_fai               // channel: [mandatory] [ fasta_fai ] => The index of the fasta reference file

    main:

    ch_versions     = Channel.empty()

    beds
        .map(
            { meta, bed, bed_gz, bed_gz_tbi ->
                [ meta, bed_gz, bed_gz_tbi ]
            }
        )
        .dump(tag: 'gzipped_beds', pretty: true)
        .set { gzipped_beds }

    //
    // Calling variants using Manta
    //

    crams
        .combine(gzipped_beds, by: 0)
        .dump(tag: 'manta_input', pretty: true)
        .set { manta_input }

    MANTA_GERMLINE(
        manta_input,
        fasta,
        fasta_fai
    )

    ch_versions = ch_versions.mix(MANTA_GERMLINE.out.versions)

    //
    // Reformat the inversions into single inverted sequence junctions
    //

    MANTA_CONVERTINVERSION(
        MANTA_GERMLINE.out.diploid_sv_vcf,
        fasta
    )

    ch_versions = ch_versions.mix(MANTA_CONVERTINVERSION.out.versions)

    MANTA_CONVERTINVERSION.out.vcf
        .combine(MANTA_CONVERTINVERSION.out.tbi, by:0)
        .map(
            { meta, vcf, tbi ->
                new_meta = meta.clone()
                new_meta.caller = "manta"
                [ new_meta, vcf, tbi ]
            }
        )
        .dump(tag: 'manta_vcfs', pretty: true)
        .set { manta_vcfs }

    emit:
    manta_vcfs
    versions = ch_versions
}
