//
// Run QDNAseq
//

include { QDNASEQ as QDNASEQ_MALE       } from '../../../modules/local/qdnaseq/main'
include { QDNASEQ as QDNASEQ_FEMALE     } from '../../../modules/local/qdnaseq/main'
include { SAMTOOLS_CONVERT              } from '../../../modules/nf-core/samtools/convert/main'
include { NGSBITS_SAMPLEGENDER          } from '../../../modules/nf-core/ngsbits/samplegender/main'

workflow BAM_VARIANT_CALLING_QDNASEQ {
    take:
        ch_crams                // channel: [mandatory] [ meta, cram, crai ] => The aligned CRAMs per sample with the regions they should be called on
        ch_fasta                // channel: [mandatory] [ meta, fasta ] => The fasta reference file
        ch_fai                  // channel: [mandatory] [ meta, fai ] => The index of the fasta reference file
        ch_qdnaseq_male         // channel: [mandatory] [ meta, qdnaseq_reference ] => The male reference to be used for qDNAseq
        ch_qdnaseq_female       // channel: [mandatory] [ meta, qdnaseq_reference ] => The female reference to be used for qDNAseq

    main:

    ch_versions     = Channel.empty()

    ch_crams
        .branch { meta, cram, crai ->
            gender: meta.gender
                return [meta, meta.gender]
            no_gender: !meta.gender
        }
        .set { ch_samplegender_input }

    NGSBITS_SAMPLEGENDER(
        ch_samplegender_input.no_gender,
        ch_fasta,
        ch_fai,
        "xy"
    )
    ch_versions = ch_versions.mix(NGSBITS_SAMPLEGENDER.out.versions.first())

    NGSBITS_SAMPLEGENDER.out.tsv
        .map { meta, tsv ->
            gender = get_gender(tsv)
            [ meta, gender ]
        }
        .mix(ch_samplegender_input.gender)
        .set { ch_genders }

    SAMTOOLS_CONVERT(
        ch_crams,
        ch_fasta.map { it[1] },
        ch_fai.map { it[1] }
    )
    ch_versions = ch_versions.mix(SAMTOOLS_CONVERT.out.versions.first())

    ch_genders
        .join(SAMTOOLS_CONVERT.out.alignment_index, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, gender, bam, bai ->
            new_meta = meta + [gender:gender]
            [ new_meta, bam, bai ]
        }
        .branch { meta, bam, bai ->
            male: meta.gender == "male"
            female: meta.gender == "female"
            other: true
        }
        .set { ch_qdnaseq_input }

    ch_qdnaseq_input.other.view { meta, bam, bai ->
        log.warn("Couldn't define the gender of sample ${meta.id}. Defaulting to male. (Specify the gender in the samplesheet to avoid this warning.)")
    }

    ch_qdnaseq_input.male
        .mix(ch_qdnaseq_input.other)
        .set { ch_males_and_others }

    QDNASEQ_MALE(
        ch_males_and_others,
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

    emit:
    qdnaseq_beds  = ch_qdnaseq_beds  // channel: [ val(meta), path(bed) ]

    versions    = ch_versions
}

def get_gender(tsv) {
    if(workflow.stubRun) {
        return "other"
    }
    split_tsv = tsv.splitCsv(sep:"\t", header:true, strip:true)
    return split_tsv[0].gender
}