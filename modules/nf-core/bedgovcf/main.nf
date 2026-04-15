nextflow.preview.types = true

process BEDGOVCF {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bedgovcf:0.1.1--h9ee0642_1'
        : 'biocontainers/bedgovcf:0.1.1--h9ee0642_1'}"

    input:
    tuple(meta: Map, bed: Path, config: Path)
    tuple(meta2: Map, fai: Path)

    output:
    vcf = tuple(meta, file("*.vcf.gz"))
    versions_bedgovcf = tuple("${task.process}", "bedgovcf", eval("bedgovcf --version 2>&1 | sed 's/^bedgovcf version //'"))
    versions_bgzip = tuple("${task.process}", "bgzip", eval('bgzip --version | head -1 | sed "s/bgzip (htslib) //"'))

    topic:
    tuple("${task.process}", "bedgovcf", eval("bedgovcf --version 2>&1 | sed 's/^bedgovcf version //'")) >> 'versions'
    tuple("${task.process}", "bgzip", eval('bgzip --version | head -1 | sed "s/bgzip (htslib) //"')) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bedgovcf \\
        ${args} \\
        --bed ${bed} \\
        --fai ${fai} \\
        --config ${config} \\
        | bgzip --stdout --threads ${task.cpus} ${args2} > ${prefix}.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    """
}
