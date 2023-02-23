//
// Merge VCFs from multiple callers
//

include { TABIX_BGZIP as UNZIP_VCFS                     } from '../../../modules/nf-core/tabix/bgzip/main'
include { TABIX_BGZIP as BGZIP_MERGED                   } from '../../../modules/nf-core/tabix/bgzip/main'
include { TABIX_TABIX                                   } from '../../../modules/nf-core/tabix/tabix/main'
include { JASMINESV                                     } from '../../../modules/nf-core/jasminesv/main'
include { BCFTOOLS_SORT                                 } from '../../../modules/nf-core/bcftools/sort/main'


workflow VCF_MERGE_JASMINE {
    take:
        vcfs                    // channel: [mandatory] [ meta, vcf ] => The gzipped called VCFs
        fasta                   // channel: [mandatory] [ fasta ] => The fasta reference file
        fasta_fai               // channel: [mandatory] [ fasta_fai ] => The index of the fasta reference file

    main:

    ch_versions     = Channel.empty()

    UNZIP_VCFS(
        vcfs
    )

    ch_versions = ch_versions.mix(UNZIP_VCFS.out.versions)

    UNZIP_VCFS.out.output
        .map { meta, vcf ->
            [ meta.findAll { !(it.key == "caller")}, vcf ]
        }
        .groupTuple(size:params.callers.tokenize(",").size)
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

    BGZIP_MERGED(
        JASMINESV.out.vcf
    )

    ch_versions = ch_versions.mix(BGZIP_MERGED.out.versions)

    BCFTOOLS_SORT(
        BGZIP_MERGED.out.output
    )

    ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions)

    TABIX_TABIX(
        BCFTOOLS_SORT.out.vcf
    )

    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

    BCFTOOLS_SORT.out.vcf
        .join(TABIX_TABIX.out.tbi)
        .set { merged_vcfs }

    emit:
    merged_vcfs
    versions = ch_versions
}
