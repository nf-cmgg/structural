//
// Annotate the VCFs
//

include { BCFTOOLS_SPLIT_BY_SVTYPE                  } from '../../../modules/local/bcftools/split_by_svtype'
include { BCFTOOLS_CONSENSUS_REHEADER               } from '../../../modules/local/bcftools/consensus_reheader'

include { ANNOTSV_ANNOTSV                           } from '../../../modules/nf-core/annotsv/annotsv/main'
include { ENSEMBLVEP_VEP                            } from '../../../modules/nf-core/ensemblvep/vep/main'
include { TABIX_TABIX as TABIX_ANNOTSV              } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_VEP                  } from '../../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_FILTER                           } from '../../../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_CONCAT                           } from '../../../modules/nf-core/bcftools/concat/main'

workflow VCF_ANNOTATE_VEP_ANNOTSV {
    take:
        ch_vcfs                                 // channel: [mandatory] [ val(meta), path(vcf), path(tbi) ] VCFs containing the called structural variants
        ch_small_variants                       // channel: [mandatory] [ val(meta), path(vcf) ] VCFs containing small variants used in AnnotSV
        ch_fasta                                // channel: [mandatory] [ val(meta), path(fasta) ] => The fasta reference file
        ch_annotsv_annotations                  // channel: [mandatory] [ val(meta), path(annotations) ] => The annotations for AnnotSV
        ch_annotsv_candidate_genes              // channel: [optional]  [ val(meta), path(candidate_genes) ]
        ch_annotsv_gene_transcripts             // channel: [optional]  [ val(meta), path(gene_transcripts) ]
        ch_vep_cache                            // channel: [optional]  [ path(cache) ] => The path to the local VEP cache
        ch_vep_extra_files                      // channel: [optional]  [ path(file1, file2, file3...) ] => The VEP extra files
        val_variant_types                       // list:    [mandatory] => The variant types
        genome                                  // string:  [mandatory] => The genome used by the variant callers
        species                                 // string:  [mandatory] => The species used by VEP
        vep_cache_version                       // integer: [mandatory] => The VEP cache version to use

    main:

    def ch_versions = Channel.empty()
    def ch_reports  = Channel.empty()

    // Run AnnotSV and VEP in parallel and merge TSV from AnnotSV with VCF from VEP during VCFanno

    BCFTOOLS_FILTER(
        ch_vcfs
    )
    ch_versions = ch_versions.mix(BCFTOOLS_FILTER.out.versions.first())

    BCFTOOLS_SPLIT_BY_SVTYPE(
        BCFTOOLS_FILTER.out.vcf.map { meta, vcf -> [ meta, vcf, [] ]}
    )
    ch_versions = ch_versions.mix(BCFTOOLS_SPLIT_BY_SVTYPE.out.versions.first())

    def ch_small_variants_types = ch_small_variants
        .combine(Channel.fromList(val_variant_types))
        .map { meta, vcf, type ->
            def new_meta = meta + [ variant_type:type ]
            [ new_meta, vcf ]
        }

    def ch_annotsv_input = BCFTOOLS_SPLIT_BY_SVTYPE.out.split_vcfs
        .join(ch_small_variants_types, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, vcfs, small_variants ->
            def new_meta = meta + [vcf_count:vcfs instanceof ArrayList ? vcfs.size() : 1]
            [ new_meta, vcfs instanceof ArrayList ? vcfs : [vcfs], small_variants ]
        }
        .transpose(by:1)
        .map { meta, vcf, small_variants ->
            def new_meta = meta + [id:vcf.name.replace(".vcf", "")]
            [ new_meta, vcf, [], small_variants ]
        }

    ANNOTSV_ANNOTSV(
        ch_annotsv_input,
        ch_annotsv_annotations,
        ch_annotsv_candidate_genes,
        [[],[]],
        ch_annotsv_gene_transcripts
    )
    ch_versions = ch_versions.mix(ANNOTSV_ANNOTSV.out.versions.first())

    def ch_consensus_reheader_input = ANNOTSV_ANNOTSV.out.vcf
        .join(ch_annotsv_input, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, annot_vcf, vcf, _tbi, _small ->
            [ meta, annot_vcf, vcf, [] ]
        }

    def val_additional_headers = [
        '##INFO=<ID=BNDrescue,Number=0,Type=Flag,Description="The other BND of this pair can be recovered"'
    ]

    BCFTOOLS_CONSENSUS_REHEADER(
        ch_consensus_reheader_input,
        [[],[]],
        val_additional_headers
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONSENSUS_REHEADER.out.versions.first())

    def ch_concat_input = BCFTOOLS_CONSENSUS_REHEADER.out.vcf
        .map { meta, vcf ->
            def new_meta = meta + [id:meta.sample] - meta.subMap("vcf_count")
            [ groupKey(new_meta, meta.vcf_count), vcf ]
        }
        .groupTuple()
        .map { meta, vcfs ->
            [ meta, vcfs, [] ]
        }

    BCFTOOLS_CONCAT(
        ch_concat_input
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions.first())

    TABIX_ANNOTSV(
        BCFTOOLS_CONCAT.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_ANNOTSV.out.versions.first())

    def ch_annotsv_output = BCFTOOLS_CONCAT.out.vcf
        .join(TABIX_ANNOTSV.out.tbi, failOnDuplicate:true, failOnMismatch:true)

    ENSEMBLVEP_VEP(
        ch_vcfs,
        genome,
        species,
        vep_cache_version,
        ch_vep_cache,
        ch_fasta,
        ch_vep_extra_files
    )
    ch_reports  = ch_reports.mix(ENSEMBLVEP_VEP.out.report)
    ch_versions = ch_versions.mix(ENSEMBLVEP_VEP.out.versions)

    TABIX_VEP(
        ENSEMBLVEP_VEP.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_VEP.out.versions)

    def ch_vep_output = ENSEMBLVEP_VEP.out.vcf
        .join(TABIX_VEP.out.tbi, failOnDuplicate:true, failOnMismatch:true)

    emit:
    vep_vcfs        = ch_vep_output     // channel: [ val(meta), path(vcf), path(tbi) ]
    annotsv_vcfs    = ch_annotsv_output // channel: [ val(meta), path(vcf), path(tbi) ]

    reports         = ch_reports
    versions        = ch_versions
}
