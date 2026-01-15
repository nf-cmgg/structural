//
// Merge VCFs from multiple callers
//

include { TABIX_TABIX               } from '../../../modules/nf-core/tabix/tabix/main'
include { JASMINESV                 } from '../../../modules/nf-core/jasminesv/main'
include { BCFTOOLS_SORT             } from '../../../modules/nf-core/bcftools/sort/main'

include { BCFTOOLS_CONSENSUS_REHEADER } from '../../../modules/local/bcftools/consensus_reheader/main'
include { FIX_CALLERS                 } from '../../../modules/local/fix_callers/main'

workflow VCF_MERGE_FAMILY_JASMINE {
    take:
        ch_vcfs     // channel: [mandatory] [ meta, vcf, tbi ] => The bgzipped called VCFs
        ch_fasta    // channel: [mandatory] [ meta, fasta ] => The fasta reference file
        ch_fai      // channel: [mandatory] [ meta, fai ] => The index of the fasta reference file

    main:

    def ch_versions     = channel.empty()

    def ch_jasmine_input = ch_vcfs
        .map { meta, vcf, tbi ->
            def new_meta = meta - meta.subMap("sample", "sex") + ["id":meta.variant_type ? "${meta.family}.${meta.variant_type}" : meta.family]
            [ groupKey(new_meta, meta.family_count), vcf, tbi ]
        }
        .groupTuple()


    JASMINESV(
        ch_jasmine_input.map { meta, vcfs, _tbis ->
            [ meta, vcfs, [], [], [] ]
        },
        ch_fasta,
        ch_fai,
        []
    )

    FIX_CALLERS(
        JASMINESV.out.vcf
    )
    ch_versions = ch_versions.mix(FIX_CALLERS.out.versions.first())

    def ch_reheader_input = FIX_CALLERS.out.vcf
        .join(ch_jasmine_input, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, merged_vcf, vcfs, tbis ->
            [ meta, merged_vcf, vcfs, tbis ]
        }
        .dump(tag:"family_reheader_input", pretty: true)

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

    def ch_vcfs_out = BCFTOOLS_SORT.out.vcf
        .join(TABIX_TABIX.out.index, failOnMismatch:true, failOnDuplicate:true)
        .map { meta, vcf, tbi ->
            def new_meta = meta + [id:meta.family]
            [ new_meta, vcf, tbi ]
        }

    emit:
    vcfs        = ch_vcfs_out    // channel: [ val(meta), path(vcf), path(tbi) ]

    versions    = ch_versions
}
