//
// Gather Sample Evidence
//
include { SVTK_STANDARDIZE } from '../../../modules/nf-core/modules/svtk/standardize/main'
include { TABIX_TABIX      } from '../../../modules/nf-core/modules/tabix/tabix/main'
include { SVTEST_VCF       } from '../../../modules/local/svtest/vcf/main'
include { SVTEST_SRFILE    } from '../../../modules/local/svtest/sr-file/main'

workflow GATHER_SAMPLE_EVIDENCE_METRICS {
    take:
        called_vcfs             // channel: [mandatory] [ meta, vcf, tbi ] => The VCFs from all used variant callers
        split_read_evidence     // channel: [mandatory] [ meta, split_read_evidence ] => The split read evidence
        paired_end_evidence     // channel: [mandatory] [ meta, paired_end_evidence ] => The paired end evidence
        fasta_fai               // channel: [mandatory] [ fasta_fai ] => The index of the fasta reference file

    main:

    ch_metrics  = Channel.empty()
    ch_versions = Channel.empty()

    SVTK_STANDARDIZE(
        called_vcfs.map({ meta, vcf, tbi -> [ meta, vcf ]}),
        fasta_fai
    )

    ch_versions = ch_versions.mix(SVTK_STANDARDIZE.out.versions)

    TABIX_TABIX(
        SVTK_STANDARDIZE.out.standardized_vcf
    )

    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

    SVTEST_VCF(
        SVTK_STANDARDIZE.out.standardized_vcf
            .combine(TABIX_TABIX.out.tbi, by:0)
            .map({ meta, vcf, tbi -> [ meta, vcf, tbi, [] ]}),
        fasta_fai
    )

    ch_metrics  = ch_metrics.mix(SVTEST_VCF.out.metrics)
    ch_versions = ch_versions.mix(SVTEST_VCF.out.versions)

    SVTEST_SRFILE(
        split_read_evidence
    )

    ch_metrics  = ch_metrics.mix(SVTEST_SRFILE.out.metrics)
    ch_versions = ch_versions.mix(SVTEST_SRFILE.out.versions)

    emit:
    metrics  = ch_metrics
    versions = ch_versions
}