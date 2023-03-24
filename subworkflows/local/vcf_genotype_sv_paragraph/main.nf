import groovy.json.JsonSlurper

//
// Merge VCFs from multiple callers
//

include { PARAGRAPH_IDXDEPTH                } from '../../../modules/nf-core/paragraph/idxdepth/main'
include { PARAGRAPH_MULTIGRMPY              } from '../../../modules/nf-core/paragraph/multigrmpy/main'
include { TABIX_TABIX as TABIX_INDIVIDUALS  } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_FAMILY       } from '../../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_MERGE                    } from '../../../modules/nf-core/bcftools/merge/main'

workflow VCF_GENOTYPE_SV_PARAGRAPH {
    take:
        vcfs                    // channel: [mandatory] [ meta, vcf, tbi ] VCFs containing the called structural variants
        crams                   // channel: [mandatory] [ meta, cram, crai ] => The CRAM files used to create the VCF files
        fasta                   // channel: [mandatory] [ fasta ] => The fasta reference file
        fai               // channel: [mandatory] [ fai ] => The index of the fasta reference file

    main:

    ch_versions     = Channel.empty()

    PARAGRAPH_IDXDEPTH(
        crams,
        fasta.map { [[], it] },
        fai.map { [[], it] }
    )
    ch_versions = ch_versions.mix(PARAGRAPH_IDXDEPTH.out.versions.first())

    PARAGRAPH_IDXDEPTH.out.depth
        .tap { meta_channel }
        .map { meta, json ->
            manifest = create_manifest(json.text, meta.id)
            [ meta, manifest ]
        }
        .collectFile() { meta, manifest ->
            [ "${meta.id}.tsv", manifest ]
        }
        .map { manifest ->
            id = manifest.baseName
            [ id, manifest ]
        }
        .join(meta_channel.map { meta, json -> [ meta.id, meta ]}, failOnMismatch:true, failOnDuplicate:true)
        .map { id, manifest, meta ->
            [ meta, manifest ]
        }
        .set { manifest }

    vcfs
        .join(crams, failOnMismatch:true, failOnDuplicate:true)
        .join(manifest, failOnMismatch:true, failOnDuplicate:true)
        .set { grmpy_input }

    PARAGRAPH_MULTIGRMPY(
        grmpy_input,
        fasta.map { [[], it] },
        fai.map { [[], it] }
    )
    ch_versions = ch_versions.mix(PARAGRAPH_MULTIGRMPY.out.versions.first())

    TABIX_INDIVIDUALS(
        PARAGRAPH_MULTIGRMPY.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_INDIVIDUALS.out.versions.first())

    PARAGRAPH_MULTIGRMPY.out.vcf
        .join(TABIX_INDIVIDUALS.out.tbi, failOnMismatch:true, failOnDuplicate:true)
        .map { meta, vcf, tbi ->
            new_meta = meta.findAll { !(it.key == "sample") } + [id:meta.family]
            [ groupKey(new_meta, meta.family_count), vcf, tbi ]
        }
        .groupTuple()
        .branch { meta, vcfs, tbis ->
            merge: vcfs.size() > 1
            dont_merge: vcfs.size() == 1
                return [ meta, vcfs ]
        }
        .set { merge_input }

    BCFTOOLS_MERGE(
        merge_input.merge,
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(BCFTOOLS_MERGE.out.versions.first())

    BCFTOOLS_MERGE.out.merged_variants
        .mix(merge_input.dont_merge)
        .set { genotyped_vcfs }

    if(!params.annotate) {
        TABIX_FAMILY(
            genotyped_vcfs
        )
        ch_versions = ch_versions.mix(TABIX_FAMILY.out.versions.first())
    }

    emit:
    genotyped_vcfs
    versions = ch_versions
}

def create_manifest(json, id) {
    if(!workflow.stubRun) {
        Map jsonMap = (Map) new JsonSlurper().parseText(json)
        initDepth = jsonMap["autosome"]["depth"]
        depth = workflow.profile.contains("test") && initDepth == 0.0 ? 0.1 : initDepth
        path = jsonMap["bam_path"]
        read_length = jsonMap["read_length"]
        sex = "unknown" //TODO: add support for sex determination
        header = "id\tpath\tdepth\tread length\tsex"
        return "${header}\n${id}\t${path}\t${depth}\t${read_length}\t${sex}".toString()
    }
    else {
        return "stub"
    }
}
