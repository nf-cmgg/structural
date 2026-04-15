nextflow.preview.types = true

process SAMTOOLS_MERGE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0'
        : 'biocontainers/samtools:1.22.1--h96c455f_0'}"

    input:
    tuple(meta: Map, input_files: List<Path>)
    tuple(meta2: Map, fasta: Path, fai: Path, gzi: Path)

    stage:
    stageAs input_files, "?/*"

    output:
    bam = tuple(meta, file("${prefix}.bam", optional: true))
    cram = tuple(meta, file("${prefix}.cram", optional: true))
    csi = tuple(meta, file("*.csi", optional: true))
    crai = tuple(meta, file("*.crai", optional: true))
    versions_samtools = tuple("${task.process}", 'samtools', eval("samtools version | sed '1!d;s/.* //'"))

    topic:
    tuple("${task.process}", 'samtools', eval("samtools version | sed '1!d;s/.* //'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def file_type = input_files.first().getExtension()
    def reference = fasta ? "--reference ${fasta}" : ""
    """
    # Note: --threads value represents *additional* CPUs to allocate (total CPUs = 1 + --threads).
    samtools \\
        merge \\
        --threads ${task.cpus - 1} \\
        ${args} \\
        ${reference} \\
        ${prefix}.${file_type} \\
        ${input_files}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    def file_type = input_files.first().getExtension()
    def index_type = file_type == "bam" ? "csi" : "crai"
    def index = args.contains("--write-index") ? "touch ${prefix}.${index_type}" : ""
    """
    touch ${prefix}.${file_type}
    ${index}
    """
}
