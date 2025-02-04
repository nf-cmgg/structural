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
        ch_svync_configs    // channel: [mandatory] [ configs ] => A list of svync configs

    main:

    def ch_versions     = Channel.empty()

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

    def ch_contigs = TABIX_BGZIPTABIX.out.gz_tbi
        .map { _meta, bed_gz, tbi ->
            [bed_gz, tbi]
        }
        .collect()

    //
    // Calling variants using Manta
    //

    def ch_manta_input = ch_crams
        .combine(ch_contigs)
        .dump(tag: 'manta_input', pretty: true)

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

    def ch_manta_svync_config = ch_svync_configs
        .map { configs ->
            configs.find { config -> config.toString().contains("manta") }
        }

    def ch_manta_vcfs = MANTA_CONVERTINVERSION.out.vcf
        .join(MANTA_CONVERTINVERSION.out.tbi, failOnDuplicate:true, failOnMismatch:true)
        
    def ch_svync_input = ch_manta_vcfs
        .combine(ch_manta_svync_config)
        .dump(tag: 'manta_vcfs', pretty: true)

    SVYNC(
        ch_svync_input
    )
    ch_versions = ch_versions.mix(SVYNC.out.versions.first())

    TABIX_TABIX(
        SVYNC.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    def ch_out_vcfs = SVYNC.out.vcf
        .join(TABIX_TABIX.out.tbi)

    emit:
    raw_vcfs    = ch_manta_vcfs // channel: [ val(meta), path(vcf), path(tbi) ]
    manta_vcfs  = ch_out_vcfs   // channel: [ val(meta), path(vcf), path(tbi) ]

    versions    = ch_versions
}
