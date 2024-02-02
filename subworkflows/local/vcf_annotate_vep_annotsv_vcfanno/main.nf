//
// Annotate the VCFs
//

include { BCFTOOLS_SPLIT_BY_SVTYPE                  } from '../../../modules/local/bcftools/split_by_svtype'
include { BCFTOOLS_CONSENSUS_REHEADER               } from '../../../modules/local/bcftools/consensus_reheader'

include { ANNOTSV_ANNOTSV                           } from '../../../modules/nf-core/annotsv/annotsv/main'
include { ENSEMBLVEP_VEP                            } from '../../../modules/nf-core/ensemblvep/vep/main'
include { VCFANNO                                   } from '../../../modules/nf-core/vcfanno/main'
include { TABIX_BGZIPTABIX as TABIX_ANNOTATED       } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_TABIX as TABIX_ANNOTSV              } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_VEP                  } from '../../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_FILTER                           } from '../../../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_FILTER as BCFTOOLS_FILTER_COMMON } from '../../../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_CONCAT                           } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_SORT                             } from '../../../modules/nf-core/bcftools/sort/main'
include { TABIX_TABIX as TABIX_FILTER               } from '../../../modules/nf-core/tabix/tabix/main'

workflow VCF_ANNOTATE_VEP_ANNOTSV_VCFANNO {
    take:
        ch_vcfs                                 // channel: [mandatory] [ val(meta), path(vcf), path(tbi) ] VCFs containing the called structural variants
        ch_small_variants                       // channel: [mandatory] [ val(meta), path(vcf) ] VCFs containing small variants used in AnnotSV
        ch_fasta                                // channel: [mandatory] [ val(meta), path(fasta) ] => The fasta reference file
        ch_fai                                  // channel: [mandatory] [ val(meta), path(fai) ] => The index of the fasta reference file
        ch_annotsv_annotations                  // channel: [mandatory] [ val(meta), path(annotations) ] => The annotations for AnnotSV
        ch_annotsv_candidate_genes              // channel: [optional]  [ val(meta), path(candidate_genes) ]
        ch_annotsv_gene_transcripts             // channel: [optional]  [ val(meta), path(gene_transcripts) ]
        ch_vep_cache                            // channel: [optional]  [ path(cache) ] => The path to the local VEP cache
        ch_vep_extra_files                      // channel: [optional]  [ path(file1, file2, file3...) ] => The VEP extra files
        ch_vcfanno_lua                          // channel: [optional]  [ path(lua) ] => The lua script to influence VCFanno
        val_vcfanno_resources                   // list:    [optional]  [ path(file1, file2, file3...) ] => The extra VCFanno files
        val_variant_types                       // list:    [mandatory] => The variant types

    main:

    ch_versions = Channel.empty()
    ch_reports  = Channel.empty()

    // Run AnnotSV and VEP in parallel and merge TSV from AnnotSV with VCF from VEP during VCFanno

    BCFTOOLS_FILTER(
        ch_vcfs.map { meta, vcf, tbi -> [ meta, vcf ]}
    )
    ch_versions = ch_versions.mix(BCFTOOLS_FILTER.out.versions.first())

    BCFTOOLS_SPLIT_BY_SVTYPE(
        BCFTOOLS_FILTER.out.vcf.map { meta, vcf -> [ meta, vcf, [] ]}
    )
    ch_versions = ch_versions.mix(BCFTOOLS_SPLIT_BY_SVTYPE.out.versions.first())

    ch_small_variants
        .combine(Channel.fromList(val_variant_types))
        .map { meta, vcf, type ->
            def new_meta = meta + [ variant_type:type ]
            [ new_meta, vcf ]
        }
        .set { ch_small_variants_types }

    BCFTOOLS_SPLIT_BY_SVTYPE.out.split_vcfs
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
        .set { ch_annotsv_input }

    ANNOTSV_ANNOTSV(
        ch_annotsv_input,
        ch_annotsv_annotations,
        ch_annotsv_candidate_genes,
        [[],[]],
        ch_annotsv_gene_transcripts
    )
    ch_versions = ch_versions.mix(ANNOTSV_ANNOTSV.out.versions.first())

    ANNOTSV_ANNOTSV.out.vcf
        .join(ch_annotsv_input, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, annot_vcf, vcf, tbi, small ->
            [ meta, annot_vcf, vcf, [] ]
        }
        .set { ch_consensus_reheader_input }

    val_additional_headers = [
        '##INFO=<ID=BNDrescue,Number=0,Type=Flag,Description="The other BND of this pair can be recovered"'
    ]

    BCFTOOLS_CONSENSUS_REHEADER(
        ch_consensus_reheader_input,
        [[],[]],
        val_additional_headers
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONSENSUS_REHEADER.out.versions.first())

    BCFTOOLS_CONSENSUS_REHEADER.out.vcf
        .map { meta, vcf ->
            def new_meta = meta + [id:meta.sample] - meta.subMap("vcf_count")
            [ groupKey(new_meta, meta.vcf_count), vcf ]
        }
        .groupTuple()
        .map { meta, vcfs ->
            [ meta, vcfs, [] ]
        }
        .set { ch_concat_input }

    BCFTOOLS_CONCAT(
        ch_concat_input
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions.first())

    TABIX_ANNOTSV(
        BCFTOOLS_CONCAT.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_ANNOTSV.out.versions.first())

    BCFTOOLS_CONCAT.out.vcf
        .join(TABIX_ANNOTSV.out.tbi, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, vcf, tbi ->
            [ meta, [vcf, tbi]]
        }
        .set { ch_annotsv_output }

    ENSEMBLVEP_VEP(
        ch_vcfs,
        params.genome,
        params.species,
        params.vep_cache_version,
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

    ENSEMBLVEP_VEP.out.vcf
        .join(TABIX_VEP.out.tbi, failOnDuplicate:true, failOnMismatch:true)
        .join(ch_annotsv_output, failOnDuplicate:true, failOnMismatch:true)
        .set { ch_vcfanno_input }

    Channel.fromList(create_vcfanno_toml(val_vcfanno_resources))
        .collectFile(name:"vcfanno.toml", newLine:true)
        .collect()
        .set { ch_vcfanno_toml }

    VCFANNO(
        ch_vcfanno_input,
        ch_vcfanno_toml,
        ch_vcfanno_lua,
        val_vcfanno_resources ? Channel.fromList(val_vcfanno_resources).collect() : []
    )
    ch_versions = ch_versions.mix(VCFANNO.out.versions)

    
    if(!params.annotations_filter) {
        TABIX_ANNOTATED(
            VCFANNO.out.vcf
        )
        ch_versions = ch_versions.mix(TABIX_ANNOTATED.out.versions)

        TABIX_ANNOTATED.out.gz_tbi
            .set { ch_annotated_vcfs }
    }

    if(params.annotations_filter) {
        BCFTOOLS_FILTER_COMMON(
            VCFANNO.out.vcf
        )
        ch_versions = ch_versions.mix(BCFTOOLS_FILTER_COMMON.out.versions.first())

        TABIX_FILTER(
            BCFTOOLS_FILTER_COMMON.out.vcf
        )
        ch_versions = ch_versions.mix(TABIX_FILTER.out.versions.first())

        ch_annotated_vcfs = BCFTOOLS_FILTER_COMMON.out.vcf.join(TABIX_FILTER.out.tbi, failOnDuplicate:true, failOnMismatch:true)
    }

    emit:
    vcfs            = ch_annotated_vcfs  // channel: [ val(meta), path(vcf), path(tbi) ]

    reports         = ch_reports
    versions        = ch_versions
}

def create_vcfanno_toml(vcfanno_resources) {
    def params_toml_files = params.vcfanno_toml ? parse_toml(params.vcfanno_toml) : [postannotation:[]]
    def assets_toml_files = parse_toml("${projectDir}/assets/vcfanno/*.toml")
    def resources = vcfanno_resources.collect { it.fileName.toString() }
    resources.add("${params.annotsv_file_name}.vcf.gz" as String)
    def output = []
    for (file_name in resources) {
        if (params_toml_files.containsKey(file_name)){
            output.add(params_toml_files[file_name])
        }
        else if (assets_toml_files.containsKey(file_name)){
            output.add(assets_toml_files[file_name])
        }
    }
    postannotation = params_toml_files.postannotation != [] ? params_toml_files.postannotation : assets_toml_files.postannotation
    if (postannotation != []){
        output.add(postannotation)
    }
    return output.flatten()
}

def parse_toml(tomls) {
    def output = [:]
    output.postannotation = []
    tomls_files = file(tomls, checkIfExists:true)
    toml_list = tomls_files instanceof LinkedList ? tomls_files : [tomls_files]
    for (toml in toml_list) {
        def info = ""
        def file = ""
        def fields = []
        for (line in toml.readLines()) {
            if (line.startsWith("#")) { continue }
            if (line == "[[annotation]]" || line == "[[postannotation]]") {
                if (info.startsWith("[[postannotation]]")) {
                    output.postannotation.add(create_toml_config(fields))
                }
                else if(info != "") {
                    output[file] = create_toml_config(fields)
                }
                if (info != "") {
                    info = ""
                    file = ""
                    fields = []
                }
                info = line
            }
            else if (line.startsWith("file")) {
                file = line.split("\"")[-1]
            }
            fields.add(line)
        }
        if (info.startsWith("[[postannotation]]")) {
            output.postannotation.add(create_toml_config(fields))
        } else {
            output[file] = create_toml_config(fields)
        }
    }
    return output
}

def create_toml_config(fields_list) {
    config = fields_list.findAll { it != "" }.join("\n")
    return "${config}\n"
}
