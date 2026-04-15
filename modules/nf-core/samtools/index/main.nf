nextflow.preview.types = true

process SAMTOOLS_INDEX {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0'
        : 'biocontainers/samtools:1.22.1--h96c455f_0'}"

    input:
    tuple(meta: Map, input: Path)

    output:
    bai = tuple(meta, file("*.bai", optional: true))
    csi = tuple(meta, file("*.csi", optional: true))
    crai = tuple(meta, file("*.crai", optional: true))
    versions_samtools = tuple("${task.process}", 'samtools', eval("samtools version | sed '1!d;s/.* //'"))

    topic:
    tuple("${task.process}", 'samtools', eval("samtools version | sed '1!d;s/.* //'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    """
    samtools \\
        index \\
        -@ ${task.cpus} \\
        ${args} \\
        ${input}
    """

    stub:
    def args = task.ext.args ?: ''
    def extension = file(input).getExtension() == 'cram'
        ? "crai"
        : args.contains("-c") ? "csi" : "bai"
    """
    touch ${input}.${extension}
    """
}
