nextflow.preview.types = true

process SMOOVE_CALL {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/smoove:0.2.8--h9ee0642_1'
        : 'biocontainers/smoove:0.2.8--h9ee0642_1'}"

    input:
    tuple(meta: Map, input: Path, index: Path, exclude_beds: Path)
    tuple(meta2: Map, fasta: Path)
    tuple(meta3: Map, fai: Path)

    output:
    vcf = tuple(meta, file("*.vcf.gz"))
    versions_smoove = tuple("${task.process}", 'smoove', eval("smoove -v |& sed -n 's/smoove version: *//p'"))

    topic:
    tuple("${task.process}", 'smoove', eval("smoove -v |& sed -n 's/smoove version: *//p'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def exclude = exclude_beds ? "--exclude ${exclude_beds}" : ""
    """
    smoove call \\
        ${args} \\
        --outdir . \\
        --name ${prefix} \\
        --fasta ${fasta} \\
        ${exclude} \\
        --processes ${task.cpus} \\
        ${input}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    """
}
