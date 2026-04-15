nextflow.preview.types = true

process SVYNC {
    tag "${input.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/svync:0.3.0--h9ee0642_0'
        : 'biocontainers/svync:0.3.0--h9ee0642_0'}"

    input:
    input: SvyncInput

    output:
    input + record(
        vcf: file("*.vcf.gz"),
        tbi: file("*.tbi"),
    )

    topic:
    tuple("${task.process}", 'svync', eval("svync --version | sed 's/svync version //'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${input.id}"

    if ("${input.vcf}" == "${prefix}.vcf.gz") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }

    """
    svync \\
        ${args} \\
        --config ${input.config} \\
        --input ${input.vcf} \\
        | bgzip --threads ${task.cpus} ${args2} > ${prefix}.vcf.gz \\
        && tabix ${args3} ${prefix}.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${input.id}"

    if ("${input.vcf}" == "${prefix}.vcf.gz") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }

    """
    echo | gzip -n > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}

record SvyncInput {
    id: String
    vcf: Path
    tbi: Path
    config: Path
}
