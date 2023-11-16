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

    main:

    ch_versions     = Channel.empty()

    SAMTOOLS_CONVERT(
        ch_crams,
        ch_fasta.map { it[1] },
        ch_fai.map { it[1] }
    )
    ch_versions = ch_versions.mix(SAMTOOLS_CONVERT.out.versions.first())

    SAMTOOLS_CONVERT.out.alignment_index
        .branch { meta, bam, bai ->
            male: meta.sex == "male"
            female: meta.sex == "female"
        }
        .set { ch_qdnaseq_input }

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

    QDNASEQ_MALE.out.bed
        .mix(QDNASEQ_FEMALE.out.bed)
        .set { ch_qdnaseq_beds }

    GAWK(
        ch_qdnaseq_beds,
        []
    )
    ch_versions = ch_versions.mix(GAWK.out.versions.first())

    GAWK.out.output
        .map { meta, bed ->
            [ meta, bed, file("${projectDir}/assets/bedgovcf/qdnaseq.yaml", checkIfExists:true)]
        }
        .set { ch_bedgovcf_input }

    BEDGOVCF(
        ch_bedgovcf_input,
        ch_fai
    )
    ch_versions = ch_versions.mix(BEDGOVCF.out.versions.first())

    TABIX_TABIX(
        BEDGOVCF.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    BEDGOVCF.out.vcf
        .join(TABIX_TABIX.out.tbi, failOnDuplicate:true, failOnMismatch:true)
        .set { ch_vcf }

    emit:
    qdnaseq_beds    = ch_qdnaseq_beds  // channel: [ val(meta), path(bed) ]
    vcf             = ch_vcf           // channel: [ val(meta), path(vcf), path(tbi) ]

    versions        = ch_versions
}

