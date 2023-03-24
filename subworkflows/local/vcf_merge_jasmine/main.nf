//
// Merge VCFs from multiple callers
//

include { TABIX_TABIX                                   } from '../../../modules/nf-core/tabix/tabix/main'
include { JASMINESV                                     } from '../../../modules/nf-core/jasminesv/main'
include { BCFTOOLS_SORT                                 } from '../../../modules/nf-core/bcftools/sort/main'

include { REHEADER_CALLED_VCFS                          } from '../../../modules/local/bcftools/reheader_called_vcfs/main'

workflow VCF_MERGE_JASMINE {
    take:
        vcfs                    // channel: [mandatory] [ meta, vcf ] => The gzipped called VCFs
        fasta                   // channel: [mandatory] [ fasta ] => The fasta reference file
        fai               // channel: [mandatory] [ fai ] => The index of the fasta reference file

    main:

    ch_versions     = Channel.empty()

    vcfs
        .map { meta, vcf ->
            [ meta.findAll { !(it.key == "caller")}, vcf ]
        }
        .groupTuple(size:params.callers.tokenize(",").size())
        .map { meta, vcfs ->
            [ meta, vcfs, [], [] ]
        }
        .dump(tag:'jasmine_input', pretty:true)
        .set { jasmine_input }

    JASMINESV(
        jasmine_input,
        [],
        [],
        []
    )

    ch_versions = ch_versions.mix(JASMINESV.out.versions)

    new_header = Channel.fromPath("${projectDir}/assets/header.txt").collect()

    REHEADER_CALLED_VCFS(
        JASMINESV.out.vcf,
        new_header,
        fai
    )
    ch_versions = ch_versions.mix(REHEADER_CALLED_VCFS.out.versions)

    BCFTOOLS_SORT(
        REHEADER_CALLED_VCFS.out.vcf
    )

    ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions)

    TABIX_TABIX(
        BCFTOOLS_SORT.out.vcf
    )

    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

    BCFTOOLS_SORT.out.vcf
        .join(TABIX_TABIX.out.tbi, failOnMismatch:true, failOnDuplicate:true)
        .set { merged_vcfs }

    emit:
    merged_vcfs
    versions = ch_versions
}
