//
// Run QDNAseq
//

include { QDNASEQ as QDNASEQ_MALE       } from '../../../modules/local/qdnaseq/main'
include { QDNASEQ as QDNASEQ_FEMALE     } from '../../../modules/local/qdnaseq/main'
include { SAMTOOLS_CONVERT              } from '../../../modules/nf-core/samtools/convert/main'
include { GAWK                          } from '../../../modules/nf-core/gawk/main'
include { BEDGOVCF                      } from '../../../modules/nf-core/bedgovcf/main'
include { TABIX_TABIX                   } from '../../../modules/nf-core/tabix/tabix/main'

workflow BAM_VARIANT_CALLING_QDNASEQ {
    take:
        ch_crams                // channel: [mandatory] [ meta, cram, crai ] => The aligned CRAMs per sample with the regions they should be called on
        ch_fasta                // channel: [mandatory] [ meta, fasta ] => The fasta reference file
        ch_fai                  // channel: [mandatory] [ meta, fai ] => The index of the fasta reference file
        ch_qdnaseq_male         // channel: [mandatory] [ meta, qdnaseq_reference ] => The male reference to be used for qDNAseq
        ch_qdnaseq_female       // channel: [mandatory] [ meta, qdnaseq_reference ] => The female reference to be used for qDNAseq
        ch_bedgovcf_configs     // channel: [mandatory] [ configs ] => A list of bedgovcf configs

    main:

    def ch_versions     = Channel.empty()

    def ch_caller_crams = ch_crams
        .map { meta, cram, crai ->
            def new_meta = meta + [caller:'qdnaseq']
            [ new_meta, cram, crai ]
        }

    SAMTOOLS_CONVERT(
        ch_caller_crams,
        ch_fasta,
        ch_fai
    )
    ch_versions = ch_versions.mix(SAMTOOLS_CONVERT.out.versions.first())

    def ch_qdnaseq_input = SAMTOOLS_CONVERT.out.bam
        .join(SAMTOOLS_CONVERT.out.bai, failOnDuplicate:true, failOnMismatch:true)
        .branch { meta, _bam, _bai ->
            male: meta.sex == "male"
            female: meta.sex == "female"
        }

    QDNASEQ_MALE(
        ch_qdnaseq_input.male,
        ch_qdnaseq_male
    )
    ch_versions = ch_versions.mix(QDNASEQ_MALE.out.versions.first())

    QDNASEQ_FEMALE(
        ch_qdnaseq_input.female,
        ch_qdnaseq_female
    )
    ch_versions = ch_versions.mix(QDNASEQ_FEMALE.out.versions.first())

    def ch_qdnaseq_beds = QDNASEQ_MALE.out.bed
        .mix(QDNASEQ_FEMALE.out.bed)

    def ch_qdnaseq_segments = QDNASEQ_MALE.out.segments
        .mix(QDNASEQ_FEMALE.out.segments)

    def ch_qdnaseq_statistics = QDNASEQ_MALE.out.statistics
        .mix(QDNASEQ_FEMALE.out.statistics)

    GAWK(
        ch_qdnaseq_beds,
        []
    )
    ch_versions = ch_versions.mix(GAWK.out.versions.first())

    def ch_qdnaseq_bedgovcf_config = ch_bedgovcf_configs
        .map { configs ->
            configs.find { config -> config.toString().contains("qdnaseq") }
        }

    def ch_bedgovcf_input = GAWK.out.output
        .combine(ch_qdnaseq_bedgovcf_config)
        .map { meta, bed, config ->
            [ meta, bed, config ]
        }

    BEDGOVCF(
        ch_bedgovcf_input,
        ch_fai
    )
    ch_versions = ch_versions.mix(BEDGOVCF.out.versions.first())

    def ch_vcf = BEDGOVCF.out.vcf
        .map { meta, vcf ->
            def new_meta = meta - meta.subMap("caller")
            [ new_meta, vcf ]
        }

    TABIX_TABIX(
        BEDGOVCF.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    emit:
    beds        = ch_qdnaseq_beds       // channel: [ val(meta), path(bed) ]
    segments    = ch_qdnaseq_segments   // channel: [ val(meta), path(bed) ]
    statistics  = ch_qdnaseq_statistics // channel: [ val(meta), path(stats) ]
    vcf         = ch_vcf                // channel: [ val(meta), path(vcf) ]

    versions    = ch_versions
}
