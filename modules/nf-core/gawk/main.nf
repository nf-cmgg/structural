nextflow.preview.types = true

process GAWK {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/gawk:5.3.0'
        : 'biocontainers/gawk:5.3.0'}"

    input:
    tuple(meta: Map, input: List<Path>)
    program_file: Path
    disable_redirect_output: Boolean

    output:
    output = tuple(meta, file("*.${suffix}"))
    versions_gawk = tuple("${task.process}", 'gawk', eval("awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//'"))

    topic:
    tuple("${task.process}", 'gawk', eval("awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    // args is used for the main arguments of the tool
    def args2 = task.ext.args2 ?: ''
    // args2 is used to specify a program when no program file has been given
    def prefix = task.ext.prefix ?: "${meta.id}"
    suffix = task.ext.suffix ?: "${input.first().getExtension()}"
    // use the first extension of the input files

    def program = program_file ? "-f ${program_file}" : "${args2}"
    def lst_gz = input.collect { file -> file.getExtension().endsWith("gz") ? file.toUriString() : null }
    def unzip = lst_gz ? "gunzip -q -f ${lst_gz.join(" ")}" : ""
    def input_cmd = input.collect { file -> file.extension == "gz" ? file.toUriString().replace(".gz", "") : file.toUriString() }.join(" ")
    def output_cmd = suffix.endsWith("gz") ? "| gzip > ${prefix}.${suffix}" : "> ${prefix}.${suffix}"
    def output = disable_redirect_output ? "" : output_cmd
    def cleanup = lst_gz ? "rm ${lst_gz.collect { file -> file.replace(".gz", "") }.join(" ")}" : ""

    input.collect { file ->
        assert file.name != "${prefix}.${suffix}" : "Input and output names are the same, set prefix in module configuration to disambiguate!"
    }

    """
    ${unzip}

    awk \\
        ${args} \\
        ${program} \\
        ${input_cmd} \\
        ${output}

    ${cleanup}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    suffix = task.ext.suffix ?: "${input.first().getExtension()}"
    def create_cmd = suffix.endsWith("gz") ? "echo '' | gzip >" : "touch"

    """
    ${create_cmd} ${prefix}.${suffix}
    """
}
