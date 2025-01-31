include { WISECONDORX_CONVERT } from '../../../modules/nf-core/wisecondorx/convert/main'
include { WISECONDORX_PREDICT } from '../../../modules/nf-core/wisecondorx/predict/main'
include { BEDGOVCF            } from '../../../modules/nf-core/bedgovcf/main'
include { TABIX_TABIX         } from '../../../modules/nf-core/tabix/tabix/main'

workflow BAM_VARIANT_CALLING_WISECONDORX {

    take:
    ch_crams            // channel: [ val(meta),  path(cram), path(crai) ]
    ch_fasta            // channel: [ val(meta2), path(fasta) ]
    ch_fai              // channel: [ val(meta3), path(fai) ]
    ch_ref              // channel: [ val(meta4), path(reference) ]
    ch_blacklist        // channel: [ val(meta5), path(blacklist) ]
    ch_bedgovcf_configs // channel: [mandatory] [ configs ] => A list of bedgovcf configs

    main:

    def ch_versions = Channel.empty()

    def ch_caller_crams = ch_crams
        .map { meta, cram, crai ->
            def new_meta = meta + [caller:'wisecondorx']
            [ new_meta, cram, crai ]
        }

    WISECONDORX_CONVERT(
        ch_caller_crams,
        ch_fasta,
        ch_fai
    )
    ch_versions = ch_versions.mix(WISECONDORX_CONVERT.out.versions.first())

    WISECONDORX_PREDICT(
        WISECONDORX_CONVERT.out.npz,
        ch_ref,
        ch_blacklist
    )
    ch_versions = ch_versions.mix(WISECONDORX_PREDICT.out.versions.first())

    def ch_wisecondorx_bedgovcf_config = ch_bedgovcf_configs
        .map { configs ->
            configs.find { config -> config.toString().contains("wisecondorx") }
        }

    def ch_bedgovcf_input = WISECONDORX_PREDICT.out.aberrations_bed
        .combine(ch_wisecondorx_bedgovcf_config)
        .map { meta, bed, config ->
            [ meta, bed, config]
        }

    BEDGOVCF(
        ch_bedgovcf_input,
        ch_fai
    )
    ch_versions = ch_versions.mix(BEDGOVCF.out.versions.first())

    def ch_vcf = BEDGOVCF.out.vcf
        .map { meta, vcf ->
            def new_meta = meta - meta.subMap("caller")
            [ new_meta, vcf ]
        }

    TABIX_TABIX(
        BEDGOVCF.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    emit:
    aberrations_bed = WISECONDORX_PREDICT.out.aberrations_bed   // channel: [ val(meta), path(bed) ]
    bins_bed        = WISECONDORX_PREDICT.out.bins_bed          // channel: [ val(meta), path(bed) ]
    segments_bed    = WISECONDORX_PREDICT.out.segments_bed      // channel: [ val(meta), path(bed) ]
    chr_statistics  = WISECONDORX_PREDICT.out.chr_statistics    // channel: [ val(meta), path(txt) ]
    chr_plots       = WISECONDORX_PREDICT.out.chr_plots         // channel: [ val(meta), [ path(png), path(png), ... ] ]
    genome_plot     = WISECONDORX_PREDICT.out.genome_plot       // channel: [ val(meta), path(png) ]
    vcf             = ch_vcf                                    // channel: [ val(meta), path(vcf) ]

    versions        = ch_versions                               // channel: path(versions.yml)
}

