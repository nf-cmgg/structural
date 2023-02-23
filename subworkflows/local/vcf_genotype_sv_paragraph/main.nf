import groovy.json.JsonSlurper

//
// Merge VCFs from multiple callers
//

include { PARAGRAPH_IDXDEPTH } from '../../../modules/nf-core/paragraph/idxdepth/main'

include { PARAGRAPH_GRMPY    } from '../../../modules/local/paragraph/grmpy/main'

workflow VCF_GENOTYPE_SV_PARAGRAPH {
    take:
        vcfs                    // channel: [mandatory] [ meta, vcf, tbi ] VCFs containing the called structural variants
        crams                   // channel: [mandatory] [ meta, cram, crai ] => The CRAM files used to create the VCF files
        fasta                   // channel: [mandatory] [ fasta ] => The fasta reference file
        fasta_fai               // channel: [mandatory] [ fasta_fai ] => The index of the fasta reference file

    main:

    ch_versions     = Channel.empty()

    PARAGRAPH_IDXDEPTH(
        crams,
        fasta.map { [[], it] },
        fasta_fai.map { [[], it] }
    )
    ch_versions = ch_versions.mix(PARAGRAPH_IDXDEPTH.out.versions)

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
        .join(meta_channel.map { meta, json -> [ meta.id, meta ]})
        .map { id, manifest, meta ->
            [ meta, manifest ]
        }
        .set { manifest }

    vcfs
        .join(crams)
        .join(manifest)
        .set { grmpy_input }

    PARAGRAPH_GRMPY(
        grmpy_input,
        fasta.map { [[], it] },
        fasta_fai.map { [[], it] }
    )
    ch_versions = ch_versions.mix(PARAGRAPH_GRMPY.out.versions)

    emit:
    genotyped_vcfs = PARAGRAPH_GRMPY.out.vcf.view()
    versions = ch_versions
}

def create_manifest(json, id) {
    Map jsonMap = (Map) new JsonSlurper().parseText(json)
    initDepth = jsonMap["autosome"]["depth"]
    depth = workflow.profile.contains("test") && initDepth == 0.0 ? 0.1 : initDepth
    path = jsonMap["bam_path"]
    read_length = jsonMap["read_length"]
    sex = "unknown" //TODO: add support for sex determination
    header = "id\tpath\tdepth\tread length\tsex"
    return "${header}\n${id}\t${path}\t${depth}\t${read_length}\t${sex}".toString()
}
