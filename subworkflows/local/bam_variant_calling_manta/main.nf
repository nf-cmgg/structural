//
// Run Manta
//
include { TABIX_BGZIPTABIX       } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { GAWK                   } from '../../../modules/nf-core/gawk/main'
include { MANTA_GERMLINE         } from '../../../modules/nf-core/manta/germline/main'
include { MANTA_CONVERTINVERSION } from '../../../modules/nf-core/manta/convertinversion/main'
include { TABIX_TABIX            } from '../../../modules/nf-core/tabix/tabix/main'
include { SVYNC                  } from '../../../modules/nf-core/svync/main'

workflow BAM_VARIANT_CALLING_MANTA {
    take:
        ch_crams            // channel: [mandatory] [ meta, cram, crai ] => The aligned CRAMs per sample with the regions they should be called on
        ch_fasta            // channel: [mandatory] [ meta, fasta ] => The fasta reference file
        ch_fai              // channel: [mandatory] [ meta, fai ] => The index of the fasta reference file
        ch_manta_config     // channel: [optional]  [ config ] => The config to pass to Manta

    main:

    ch_versions     = Channel.empty()

    //
    // Create a contigs BED file
    //

    GAWK(
        ch_fai,
        []
    )
    ch_versions = ch_versions.mix(GAWK.out.versions)

    TABIX_BGZIPTABIX(
        GAWK.out.output
    )
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions)

    TABIX_BGZIPTABIX.out.gz_tbi
        .map { meta, bed_gz, tbi ->
            [bed_gz, tbi]
        }
        .collect()
        .set { ch_contigs }

    //
    // Calling variants using Manta
    //

    ch_crams
        .combine(ch_contigs)
        .dump(tag: 'manta_input', pretty: true)
        .set { ch_manta_input }

    MANTA_GERMLINE(
        ch_manta_input,
        ch_fasta,
        ch_fai,
        ch_manta_config
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
        .map{ meta, vcf, tbi ->
            new_meta = meta + [caller:"manta"]
            [ new_meta, vcf, tbi, file("${projectDir}/assets/svync/manta.yaml") ]
        }
        .dump(tag: 'manta_vcfs', pretty: true)
        .set { ch_manta_vcfs }

    SVYNC(
        ch_manta_vcfs
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
    manta_vcfs  = ch_out_vcfs  // channel: [ val(meta), path(vcf), path(tbi) ]

    versions    = ch_versions
}
