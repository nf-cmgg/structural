//
// Run Gridss
//

include { GRIDSS_GRIDSS             } from '../../../modules/nf-core/gridss/gridss/main'
include { ESTIMATE_READ_LENGTH      } from '../../../modules/local/estimate_read_length/main'
include { TABIX_TABIX               } from '../../../modules/nf-core/tabix/tabix/main'
include { VIOLA                     } from '../../../modules/local/viola/main'


workflow BAM_VARIANT_CALLING_GRIDSS {
    take:
        ch_crams       // channel: [mandatory] [ meta, cram, crai ] => The aligned CRAMs per sample with the regions they should be called on
        ch_fasta       // channel: [mandatory] [ meta, fasta ] => The fasta reference file
        ch_fai         // channel: [mandatory] [ meta, fai ] => The index of the fasta reference file
        ch_bwa_index   // channel: [mandatory] [ meta, index ] => The BWA MEM index

    main:

    ch_versions     = Channel.empty()

    ESTIMATE_READ_LENGTH(
        ch_crams,
        ch_fasta,
        ch_fai
    )
    ch_versions = ch_versions.mix(ESTIMATE_READ_LENGTH.out.versions.first())

    GRIDSS_GRIDSS(
        ch_crams.map {meta, cram, crai -> [meta, cram, []]},
        ch_fasta,
        ch_fai,
        ch_bwa_index
    )
    ch_versions = ch_versions.mix(GRIDSS_GRIDSS.out.versions.first())

    GRIDSS_GRIDSS.out.vcf
        .join(ESTIMATE_READ_LENGTH.out.read_length, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, vcf, read_length ->
            new_meta = meta + [read_length:read_length]
            [ new_meta, vcf ]
        }
        .set { ch_viola_input }

    VIOLA(
        ch_viola_input,
        "gridss"
    )
    ch_versions = ch_versions.mix(VIOLA.out.versions.first())

    TABIX_TABIX(
        VIOLA.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    VIOLA.out.vcf
        .join(TABIX_TABIX.out.tbi, failOnMismatch:true, failOnDuplicate:true)
        .map{ meta, vcf, tbi ->
            new_meta = (meta - meta.subMap("read_length")) + [caller:"gridss"]
            [ new_meta, vcf, tbi ]
        }
        .dump(tag: 'gridss_vcfs', pretty: true)
        .set { ch_gridss_vcfs }

    emit:
    gridss_vcfs = ch_gridss_vcfs // channel: [ val(meta), path(vcf), path(tbi) ]
    versions    = ch_versions
}
