nextflow.preview.types = true

process SAMTOOLS_CONVERT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0'
        : 'biocontainers/samtools:1.22.1--h96c455f_0'}"

    input:
    tuple(meta: Map, input: Path, index: Path)
    tuple(meta2: Map, fasta: Path, fai: Path)

    output:
    bam = tuple(meta, file("*.bam", optional: true))
    cram = tuple(meta, file("*.cram", optional: true))
    bai = tuple(meta, file("*.bai", optional: true))
    crai = tuple(meta, file("*.crai", optional: true))
    versions_samtools = tuple("${task.process}", 'samtools', eval("samtools version | sed '1!d;s/.* //'"))

    topic:
    tuple("${task.process}", 'samtools', eval("samtools version | sed '1!d;s/.* //'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_extension = input.getExtension() == "bam" ? "cram" : "bam"

    """
    samtools view \\
        --threads ${task.cpus} \\
        --reference ${fasta} \\
        ${args} \\
        ${input} \\
        -o ${prefix}.${output_extension}

    samtools index -@${task.cpus} ${prefix}.${output_extension}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_extension = input.getExtension() == "bam" ? "cram" : "bam"
    def index_extension = input.getExtension() == "bam" ? "crai" : "bai"

    """
    touch ${prefix}.${output_extension}
    touch ${prefix}.${output_extension}.${index_extension}
    """
}
