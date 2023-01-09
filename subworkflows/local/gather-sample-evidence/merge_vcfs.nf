//
// Merge VCFs from multiple callers
//

include { TABIX_BGZIP                                   } from '../../../modules/nf-core/tabix/bgzip/main'
include { JASMINESV                                     } from '../../../modules/nf-core/jasminesv/main'

workflow MERGE_VCFS {
    take:
        vcfs                    // channel: [mandatory] [ meta, vcf ] => The gzipped called VCFs
        fasta                   // channel: [mandatory] [ fasta ] => The fasta reference file
        fasta_fai               // channel: [mandatory] [ fasta_fai ] => The index of the fasta reference file

    main:

    ch_versions     = Channel.empty()

    TABIX_BGZIP(
        vcfs
    )

    TABIX_BGZIP.out.output
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

    emit:
    merged_vcfs = JASMINESV.out.vcf
    versions = ch_versions
}
