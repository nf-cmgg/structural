nextflow.preview.types = true

process SVTOOLS_VCFTOBEDPE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/svtools:0.5.1--py_0'
        : 'biocontainers/svtools:0.5.1--py_0'}"

    input:
    tuple(meta: Map, vcf: Path)

    output:
    bedpe = tuple(meta, file("*.bedpe"))
    versions_svtools = tuple("${task.process}", 'svtools', eval("svtools --version |& sed 's/svtools //'"))

    topic:
    tuple("${task.process}", 'svtools', eval("svtools --version |& sed 's/svtools //'")) >> 'versions'

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    svtools vcftobedpe \\
        ${args} \\
        --input ${vcf} \\
        --output ${prefix}.bedpe \\
        --tempdir ./tmp
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bedpe
    """
}
