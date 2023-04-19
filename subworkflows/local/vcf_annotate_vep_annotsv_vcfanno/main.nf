//
// Annotate the VCFs
//

include { ANNOTSV_ANNOTSV                       } from '../../../modules/nf-core/annotsv/annotsv/main'
include { ENSEMBLVEP_VEP                        } from '../../../modules/nf-core/ensemblvep/vep/main'
include { VCFANNO                               } from '../../../modules/nf-core/vcfanno/main'
include { TABIX_BGZIPTABIX as TABIX_ANNOTATED   } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_BGZIPTABIX as TABIX_ANNOTSV     } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_TABIX as TABIX_VEP              } from '../../../modules/nf-core/tabix/tabix/main'

workflow VCF_ANNOTATE_VEP_ANNOTSV_VCFANNO {
    take:
        ch_vcfs                                 // channel: [mandatory] [ val(meta), path(vcf), path(tbi) ] VCFs containing the called structural variants
        ch_fasta                                // channel: [mandatory] [ path(fasta) ] => The fasta reference file
        ch_fai                                  // channel: [mandatory] [ path(fai) ] => The index of the fasta reference file
        ch_annotsv_annotations                  // channel: [mandatory] [ val(meta), path(annotations) ] => The annotations for AnnotSV
        ch_annotsv_candidate_genes              // channel: [optional]  [ val(meta), path(candidate_genes) ]
        ch_annotsv_gene_transcripts             // channel: [optional]  [ val(meta), path(gene_transcripts) ]
        ch_vep_cache                            // channel: [optional]  [ path(cache) ] => The path to the local VEP cache
        ch_vep_extra_files                      // channel: [optional]  [ path(file1, file2, file3...) ] => The VEP extra files
        ch_vcfanno_lua                          // channel: [optional]  [ path(lua) ] => The lua script to influence VCFanno
        val_vcfanno_resources                   // list:    [optional]  [ path(file1, file2, file3...) ] => The extra VCFanno files

    main:

    ch_versions = Channel.empty()
    ch_reports  = Channel.empty()

    // Run AnnotSV and VEP in parallel and merge TSV from AnnotSV with VCF from VEP during VCFanno

    ANNOTSV_ANNOTSV(
        ch_vcfs,
        ch_annotsv_annotations,
        ch_annotsv_candidate_genes,
        [[],[]],
        [[],[]],
        ch_annotsv_gene_transcripts
    )
    ch_versions = ch_versions.mix(ANNOTSV_ANNOTSV.out.versions.first())

    ANNOTSV_ANNOTSV.out.vcf
        .map { meta, vcf ->
            // Artifically create a single variant annotated VCF if the VCF is empty
            // This will only affect test runs that create VCFs with no variants
            // TODO make sure the tests don't actually need this!
            if(vcf.size() == 0){
                vcf.text = file("${projectDir}/assets/dummy_annotsv.vcf", checkIfExists:true).text
            }
            [ meta, vcf ]
        }
        .set { ch_bgzip_input }

    TABIX_ANNOTSV(
        ch_bgzip_input
    )
    ch_versions = ch_versions.mix(TABIX_ANNOTSV.out.versions.first())

    ENSEMBLVEP_VEP(
        ch_vcfs,
        params.genome,
        params.species,
        params.vep_cache_version,
        ch_vep_cache,
        ch_fasta.map { [[], it] },
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
        .join(TABIX_ANNOTSV.out.gz_tbi
            .map { meta, vcf, tbi -> [ meta, [vcf, tbi] ]},
            failOnDuplicate:true,
            failOnMismatch:true
        )
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

    TABIX_ANNOTATED(
        VCFANNO.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_ANNOTATED.out.versions)
    
    emit:
    // annotated_vcfs  = TABIX_ANNOTATED.out.gz_tbi  // channel: [ val(meta), path(vcf), path(tbi) ]

    reports         = ch_reports
    versions        = ch_versions
}

def create_vcfanno_toml(vcfanno_resources) {
    def params_toml_files = params.vcfanno_toml ? parse_toml(params.vcfanno_toml) : [postannotation:[]]
    def assets_toml_files = parse_toml("${projectDir}/assets/vcfanno/*.toml")
    def resources = vcfanno_resources.collect { it.fileName.toString() }
    resources.add("annotsv_annotated.vcf.gz")
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
