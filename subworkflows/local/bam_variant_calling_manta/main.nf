//
// Run Manta
//
include { MANTA_GERMLINE         } from '../../../modules/nf-core/manta/germline/main'
include { MANTA_CONVERTINVERSION } from '../../../modules/nf-core/manta/convertinversion/main'
include { BCFTOOLS_REHEADER      } from '../../../modules/nf-core/bcftools/reheader/main'

workflow BAM_VARIANT_CALLING_MANTA {
    take:
        ch_crams    // channel: [mandatory] [ meta, cram, crai ] => The aligned CRAMs per sample with the regions they should be called on
        ch_beds     // channel: [optional]  [ meta, bed, bed_gz, bed_gz_tbi ] => A channel containing the normal BED, the bgzipped BED and its index file
        ch_fasta    // channel: [mandatory] [ fasta ] => The fasta reference file
        ch_fai      // channel: [mandatory] [ fai ] => The index of the fasta reference file

    main:

    ch_versions     = Channel.empty()

    ch_beds
        .map(
            { meta, bed, bed_gz, bed_gz_tbi ->
                [ meta, bed_gz, bed_gz_tbi ]
            }
        )
        .dump(tag: 'gzipped_beds', pretty: true)
        .set { ch_gzipped_beds }

    //
    // Calling variants using Manta
    //

    ch_crams
        .join(ch_gzipped_beds, failOnDuplicate:true, failOnMismatch:true)
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
        ch_fasta
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
