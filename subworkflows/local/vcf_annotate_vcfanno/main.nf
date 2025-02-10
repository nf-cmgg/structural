//
// Annotate the VCFs
//

include { VCFANNO           } from '../../../modules/nf-core/vcfanno/main'
include { BCFTOOLS_FILTER   } from '../../../modules/nf-core/bcftools/filter/main'

workflow VCF_ANNOTATE_VCFANNO {
    take:
        ch_vcfs                                 // channel: [mandatory] [ val(meta), path(vcf), path(tbi) ] VCFs containing the called structural variants
        ch_sample_specific_resources            // channel: [optional]  [ val(meta), path(vcf), path(tbi) ] Files containing resources that are sample-specific
        ch_vcfanno_lua                          // channel: [optional]  [ path(lua) ] => The lua script to influence VCFanno
        val_vcfanno_resources                   // list:    [optional]  [ path(file1, file2, file3...) ] => The extra VCFanno files
        filter                                  // string:  [optional]  => A filter pattern to use after annotating
        vcfanno_toml                            // file:    [optional]  => A vcfanno TOML config
        default_vcfanno_tomls                   // list:    [mandatory] => A list of default vcfanno configs to be concatenated with the input TOML
        annotate                                // boolean: [mandatory] => Whether or not to run the full annotation or only the specified annotations

    main:

    def ch_versions = Channel.empty()
    def ch_reports  = Channel.empty()

    def ch_vcfanno_toml = Channel.empty()
    def ch_vcfanno_input = Channel.empty()
    if(annotate) {
        ch_vcfanno_toml = Channel.fromList(create_vcfanno_toml(val_vcfanno_resources, vcfanno_toml, default_vcfanno_tomls))
            .collectFile(name:"vcfanno.toml", newLine:true)
            .collect()
        def ch_collected_specific_resources = ch_sample_specific_resources.map { entry ->
            [ entry[0], entry[1..-1] ]
        }
        ch_vcfanno_input = ch_vcfs.join(ch_collected_specific_resources, failOnMismatch:true, failOnDuplicate:true)
    } else {
        ch_vcfanno_toml = Channel.fromPath(vcfanno_toml)
        ch_vcfanno_input = ch_vcfs.map { meta, vcf, tbi ->
            [ meta, vcf, tbi, [] ]
        }
    }

    VCFANNO(
        ch_vcfanno_input,
        ch_vcfanno_toml,
        ch_vcfanno_lua,
        val_vcfanno_resources ? Channel.fromList(val_vcfanno_resources).collect() : []
    )
    ch_versions = ch_versions.mix(VCFANNO.out.versions)

    def ch_vcfanno_output = VCFANNO.out.vcf.join(VCFANNO.out.tbi, failOnDuplicate:true, failOnMismatch:true)

    def ch_annotated_vcfs = Channel.empty()
    if(filter) {
        BCFTOOLS_FILTER(
            VCFANNO.out.vcf
        )
        ch_versions = ch_versions.mix(BCFTOOLS_FILTER.out.versions.first())

        ch_annotated_vcfs = BCFTOOLS_FILTER.out.vcf.join(BCFTOOLS_FILTER.out.tbi, failOnDuplicate:true, failOnMismatch:true)
    } else {
        ch_annotated_vcfs = ch_vcfanno_output
    }

    emit:
    vcfs            = ch_annotated_vcfs  // channel: [ val(meta), path(vcf), path(tbi) ]

    reports         = ch_reports
    versions        = ch_versions
}

def create_vcfanno_toml(vcfanno_resources, input_vcfanno_toml, List<Path> vcfanno_defaults) {
    def vcfanno_toml = input_vcfanno_toml ? parse_toml(input_vcfanno_toml) : [:]
    def default_tomls = parse_toml(vcfanno_defaults)
    def resources = vcfanno_resources.collect { Path resource -> resource.fileName.toString() }
    resources.add("annotsv_annotated.vcf.gz" as String)
    def output = []
    resources.each { file_name ->
        if (vcfanno_toml.containsKey(file_name)){
            output.add(vcfanno_toml[file_name])
        }
        else if (default_tomls.containsKey(file_name)){
            output.add(default_tomls[file_name])
        }
    }
    def postannotation = vcfanno_toml.postannotation != [] ? vcfanno_toml.postannotation : default_tomls.postannotation
    if (postannotation != []){
        output.add(postannotation)
    }
    return output.flatten()
}

def parse_toml(tomls) {
    def output = [:]
    output.postannotation = []
    def toml_list = tomls instanceof List ? tomls : [tomls]
    toml_list.each { toml ->
        def info = ""
        def file = ""
        def fields = []
        toml.readLines().each { line ->
            if (line.startsWith("#")) { return }
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
    def config = fields_list.findAll { field -> field != "" }.join("\n")
    return "${config}\n"
}
