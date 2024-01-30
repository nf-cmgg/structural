//
// Merge VCFs from multiple callers
//

include { TABIX_TABIX               } from '../../../modules/nf-core/tabix/tabix/main'
include { JASMINESV                 } from '../../../modules/nf-core/jasminesv/main'
include { BCFTOOLS_SORT             } from '../../../modules/nf-core/bcftools/sort/main'

include { BCFTOOLS_CONSENSUS_REHEADER } from '../../../modules/local/bcftools/consensus_reheader/main'

workflow VCF_MERGE_FAMILY_JASMINE {
    take:
        ch_vcfs     // channel: [mandatory] [ meta, vcf ] => The bgzipped called VCFs
        ch_fasta    // channel: [mandatory] [ meta, fasta ] => The fasta reference file
        ch_fai      // channel: [mandatory] [ meta, fai ] => The index of the fasta reference file

    main:

    ch_versions     = Channel.empty()

    ch_vcfs
        .filter { it[0].family_count > 1 }
        .map { meta, vcf, tbi ->
            def new_meta = meta - meta.subMap("sample", "sex") + ["id":meta.variant_type ? "${meta.family}.${meta.variant_type}" : meta.family]
            [ groupKey(new_meta, meta.family_count), vcf, tbi ]
        }
        .groupTuple()
        .tap { ch_consensus_reheader_input }
        .map { meta, vcfs, tbis -> 
            [ meta.id, meta, vcfs ]
        }
        .tap { ch_meta_file_list }
        .map { id, meta, vcfs ->
            [ "${id}_list.txt", vcfs.collect { it.baseName }.join("\n") ]
        }
        .collectFile()
        .map { 
            def id = it.name.replaceAll("_list.txt\$", "")
            [ id, it ]
        }
        .join(ch_meta_file_list, failOnMismatch:true, failOnDuplicate:true)
        .map { id, file_list, meta, vcfs ->
            [ meta, vcfs, [], [], file_list ]
        }
        .set { ch_jasmine_input }

    JASMINESV(
        ch_jasmine_input,
        ch_fasta.map{it[1]},
        ch_fai.map{it[1]},
        []
    )
    ch_versions = ch_versions.mix(JASMINESV.out.versions.first())

    JASMINESV.out.vcf
        .join(ch_consensus_reheader_input, failOnDuplicate:true, failOnMismatch:true)
        .dump(tag:"vcf_merge_family_jasmine,bcftools_consensus_reheader", pretty: true)
        .set { ch_reheader_input }

    BCFTOOLS_CONSENSUS_REHEADER(
        ch_reheader_input,
        ch_fai
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
