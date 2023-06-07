//
// Merge VCFs from multiple callers
//

include { TABIX_TABIX                                   } from '../../../modules/nf-core/tabix/tabix/main'
include { JASMINESV                                     } from '../../../modules/nf-core/jasminesv/main'
include { BCFTOOLS_SORT                                 } from '../../../modules/nf-core/bcftools/sort/main'

include { REHEADER_CALLED_VCFS                          } from '../../../modules/local/bcftools/reheader_called_vcfs/main'

workflow VCF_MERGE_JASMINE {
    take:
        ch_vcfs     // channel: [mandatory] [ meta, vcf ] => The bgzipped called VCFs
        ch_fasta    // channel: [mandatory] [ meta, fasta ] => The fasta reference file
        ch_fai      // channel: [mandatory] [ meta, fai ] => The index of the fasta reference file

    main:

    ch_versions     = Channel.empty()

    ch_vcfs
        .map { meta, vcf ->
            new_meta = meta - meta.subMap("caller")
            [ new_meta, vcf ]
        }
        .groupTuple(size:params.callers.tokenize(",").size())
        .map { meta, vcfs ->
            [ meta, vcfs, [], [] ]
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
        .collect()
        .set { ch_new_header }

    REHEADER_CALLED_VCFS(
        JASMINESV.out.vcf,
        ch_new_header,
        ch_fai
    )
    ch_versions = ch_versions.mix(REHEADER_CALLED_VCFS.out.versions.first())

    BCFTOOLS_SORT(
        REHEADER_CALLED_VCFS.out.vcf
    )

    ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions.first())

    TABIX_TABIX(
        BCFTOOLS_SORT.out.vcf
    )

    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    BCFTOOLS_SORT.out.vcf
        .join(TABIX_TABIX.out.tbi, failOnMismatch:true, failOnDuplicate:true)
        .set { ch_merged_vcfs }

    emit:
    merged_vcfs = ch_merged_vcfs    // channel: [ val(meta), path(vcf), path(tbi) ]

    versions    = ch_versions
}
