//
// Merge VCFs from multiple callers
//

include { TABIX_TABIX               } from '../../../modules/nf-core/tabix/tabix/main'
include { JASMINESV                 } from '../../../modules/nf-core/jasminesv/main'
include { BCFTOOLS_SORT             } from '../../../modules/nf-core/bcftools/sort/main'

include { BCFTOOLS_CONSENSUS_REHEADER } from '../../../modules/local/bcftools/consensus_reheader/main'
include { FIX_CALLERS                 } from '../../../modules/local/fix_callers/main'

workflow VCF_MERGE_CALLERS_JASMINE {
    take:
        ch_vcfs     // channel: [mandatory] [ meta, vcf, tbi ] => The bgzipped called VCFs
        ch_fasta    // channel: [mandatory] [ meta, fasta ] => The fasta reference file
        ch_fai      // channel: [mandatory] [ meta, fai ] => The index of the fasta reference file
        val_callers // value:   [mandatory] => The callers used
        val_type    // value:   [mandatory] => The type of variants

    main:

    ch_versions     = Channel.empty()

    ch_vcfs
        .map { meta, vcf, tbi ->
            new_meta = meta - meta.subMap("caller") + ["variant_type":val_type]
            [ new_meta, vcf, tbi ]
        }
        .groupTuple(size:val_callers.size())
        .tap { ch_consensus_reheader_input }
        .map { meta, vcfs, tbis ->
            [ meta, vcfs, [], [], [] ]
        }
        .dump(tag:'jasmine_input', pretty:true)
        .set { ch_jasmine_input }

    JASMINESV(
        ch_jasmine_input,
        ch_fasta.map{it[1]},
        ch_fai.map{it[1]},
        []
    )
    ch_versions = ch_versions.mix(JASMINESV.out.versions.first())

    FIX_CALLERS(
        JASMINESV.out.vcf
    )
    ch_versions = ch_versions.mix(FIX_CALLERS.out.versions.first())

    FIX_CALLERS.out.vcf
        .join(ch_consensus_reheader_input, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, vcf, vcfs, tbis ->
            [ meta, vcf, vcfs, [] ]
        }
        .dump(tag:"caller_reheader_input", pretty: true)
        .set { ch_reheader_input }

    BCFTOOLS_CONSENSUS_REHEADER(
        ch_reheader_input,
        ch_fai,
        []
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONSENSUS_REHEADER.out.versions.first())

    BCFTOOLS_SORT(
        BCFTOOLS_CONSENSUS_REHEADER.out.vcf
    )
    ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions.first())

    TABIX_TABIX(
        BCFTOOLS_SORT.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    BCFTOOLS_SORT.out.vcf
        .join(TABIX_TABIX.out.tbi, failOnMismatch:true, failOnDuplicate:true)
        .set { ch_vcfs_out }

    emit:
    vcfs        = ch_vcfs_out    // channel: [ val(meta), path(vcf), path(tbi) ]

    versions    = ch_versions
}
