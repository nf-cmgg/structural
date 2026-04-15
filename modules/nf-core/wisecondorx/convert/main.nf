nextflow.preview.types = true

process WISECONDORX_CONVERT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/wisecondorx:1.2.9--pyhdfd78af_0'
        : 'biocontainers/wisecondorx:1.2.9--pyhdfd78af_0'}"

    input:
    tuple(meta: Map, bam: Path, bai: Path)
    tuple(meta2: Map, fasta: Path)
    tuple(meta3: Map, fasta_fai: Path)

    output:
    npz = tuple(meta, file("*.npz"))
    versions_wisecondorx = tuple("${task.process}", 'wisecondorx', eval("pip list |& sed -n 's/wisecondorx *//p'"))

    topic:
    tuple("${task.process}", 'wisecondorx', eval("pip list |& sed -n 's/wisecondorx *//p'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--reference ${fasta}" : ""

    """
    WisecondorX convert \\
        ${bam} \\
        ${prefix}.npz \\
        ${reference} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.npz
    """
}
