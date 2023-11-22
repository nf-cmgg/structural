//
// Merge VCFs from multiple callers
//

include { TABIX_TABIX                                   } from '../../../modules/nf-core/tabix/tabix/main'
include { JASMINESV                                     } from '../../../modules/nf-core/jasminesv/main'
include { BCFTOOLS_SORT                                 } from '../../../modules/nf-core/bcftools/sort/main'

include { BCFTOOLS_REHEADER                          } from '../../../modules/nf-core/bcftools/reheader/main'

workflow VCF_MERGE_CALLERS_JASMINE {
    take:
        ch_vcfs     // channel: [mandatory] [ meta, vcf ] => The bgzipped called VCFs
        ch_fasta    // channel: [mandatory] [ meta, fasta ] => The fasta reference file
        ch_fai      // channel: [mandatory] [ meta, fai ] => The index of the fasta reference file
        val_callers // value:   [mandatory] => The callers used
        val_type    // value:   [mandatory] => The type of variants

    main:

    ch_versions     = Channel.empty()

    if(val_callers.size() > 1){
        ch_vcfs
            .map { meta, vcf ->
                new_meta = meta - meta.subMap("caller") + ["variant_type":val_type]
                [ new_meta, vcf ]
            }
            .groupTuple(size:val_callers.size())
            .map { meta, vcfs ->
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

        Channel.fromPath("${projectDir}/assets/header.txt")
            .map { meta, header ->
                [ meta, header, [] ]
            }
            .collect()
            .set { ch_new_header }

        BCFTOOLS_REHEADER(
            JASMINESV.out.vcf.combine(ch_new_header),
            ch_fai
        )
        ch_versions = ch_versions.mix(BCFTOOLS_REHEADER.out.versions.first())

        BCFTOOLS_SORT(
            BCFTOOLS_REHEADER.out.vcf
        )
        ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions.first())

        BCFTOOLS_SORT.out.vcf
            .set { ch_merged_vcfs }
    } else {
        ch_vcfs
            .map { meta, vcf ->
                new_meta = meta - meta.subMap("caller") + ["variant_type":val_type]
                [ new_meta, vcf ]
            }
            .set { ch_merged_vcfs }
    }

    TABIX_TABIX(
        ch_merged_vcfs
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    ch_merged_vcfs
        .join(TABIX_TABIX.out.tbi, failOnMismatch:true, failOnDuplicate:true)
        .set { ch_vcfs_out }

    emit:
    vcfs        = ch_vcfs_out    // channel: [ val(meta), path(vcf), path(tbi) ]

    versions    = ch_versions
}
