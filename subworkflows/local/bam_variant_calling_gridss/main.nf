//
// Run Gridss
//

include { GRIDSS_GRIDSS             } from '../../../modules/nf-core/gridss/gridss/main'
include { TABIX_TABIX               } from '../../../modules/nf-core/tabix/tabix/main'
include { VIOLA                     } from '../../../modules/local/viola/main'
include { BCFTOOLS_SORT             } from '../../../modules/nf-core/bcftools/sort/main'


workflow BAM_VARIANT_CALLING_GRIDSS {
    take:
        ch_crams       // channel: [mandatory] [ meta, cram, crai ] => The aligned CRAMs per sample with the regions they should be called on
        ch_fasta       // channel: [mandatory] [ meta, fasta ] => The fasta reference file
        ch_fai         // channel: [mandatory] [ meta, fai ] => The index of the fasta reference file
        ch_bwa_index   // channel: [mandatory] [ meta, index ] => The BWA MEM index

    main:

    ch_versions     = Channel.empty()

    GRIDSS_GRIDSS(
        ch_crams.map {meta, cram, crai -> [meta, cram, []]},
        ch_fasta,
        ch_fai,
        ch_bwa_index
    )
    ch_versions = ch_versions.mix(GRIDSS_GRIDSS.out.versions.first())

    VIOLA(
        GRIDSS_GRIDSS.out.vcf,
        "gridss"
    )
    ch_versions = ch_versions.mix(VIOLA.out.versions.first())

    BCFTOOLS_SORT(
        VIOLA.out.vcf
    )
    ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions.first())

    TABIX_TABIX(
        BCFTOOLS_SORT.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    BCFTOOLS_SORT.out.vcf
        .join(TABIX_TABIX.out.tbi, failOnMismatch:true, failOnDuplicate:true)
        .map{ meta, vcf, tbi ->
            new_meta = meta + [caller:"gridss"]
            [ new_meta, vcf, tbi ]
        }
        .dump(tag: 'gridss_vcfs', pretty: true)
        .set { ch_gridss_vcfs }

    emit:
    gridss_vcfs = ch_gridss_vcfs // channel: [ val(meta), path(vcf), path(tbi) ]
    versions    = ch_versions
}
